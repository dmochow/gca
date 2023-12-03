import numpy as np
import scipy

def gc_block_covariance(X,L):
    """
    generate the spatiotemporal covariance matrices needed to run GCA
    X is a samples-by-sensors observation array
    L is a scalar denoting the max lag in the VAR model

    returns blkR1, Rtilde (see paper)
    """

    n_samples, n_dims = X.shape

    ## Rtilde

    # initialize block covariance
    Rtilde = np.zeros(shape=(L*n_dims,L*n_dims))

    for i in range(L):
        for j in range(L):
            Xl = np.roll(X, i, axis=0)
            Xr = np.roll(X, j, axis=0)
            Xl[0:i,:]=0
            Xr[0:j,:]=0

            Rij = (1/n_samples) * (Xl.T@Xr) # biased
            #Rij = ( 1 / (n_samples-np.max([i,j])) ) * (Xl.T@Xr) # unbiased

            Rtilde[i*n_dims:(i+1)*n_dims,j*n_dims:(j+1)*n_dims]=Rij

    ## blkR1

    # lag 1
    Xp = np.roll(X,1,axis=0)
    Xp[0:1,:]=0
    blkR1 = (1/n_samples)*(X.T@Xp)
    #blkR1 = (1/(n_samples-1))*(X.T@Xp)

    # add remaining lags (block_diag adds a row for an empty matrix)
    for i in range(2,L+1):
        Xp = np.roll(X,i,axis=0)
        Xp[0:i,:]=0
        blkR1 = scipy.linalg.block_diag(blkR1, (1/n_samples)*(X.T@Xp))
        #blkR1 = scipy.linalg.block_diag(blkR1, (1/(n_samples-i))*(X.T@Xp))

    return blkR1, Rtilde


def phi_r(v,blkR1,Rtilde,L):
    """
    closed-form representation of the mmse of the "reduced" model in the GCA objective

    :param v: spatial filter
    :param blkR1: block covariance 1 see paper
    :param Rtilde: block covariance 2 see paper
    :param L: max time lag in G-causality
    :return: mmse of reduced model \Phi_r in paper
    """

    # get dimensionality (num of channels)
    D=v.shape[0]

    # extract the lag 0 covariance
    R0 = Rtilde[0:D,0:D]

    kronILv = np.kron(np.eye(L),v)
    kron1Lv = np.kron(np.ones((L,1)) , v  )

    q = kronILv.T @ blkR1 @ kron1Lv
    Q = kronILv.T @ Rtilde @ kronILv

    res = v.T @ R0 @ v - q.T @ np.linalg.pinv(Q) @ q

    return res

def phi_f(v,w,blkR1,Rtilde,L):
    """
    closed-form representation of the mmse of the "full" model in the GCA objective

    :param w: spatial filter of driving signal
    :param v: spatial filter of driven signal
    :param blkR1: block covariance 1 see paper
    :param Rtilde: block covariance 2 see paper
    :param L: max time lag in G-causality
    :return: mmse of full model \Phi_f in paper
    """

    # get dimensionality (num of channels)
    D=v.shape[0]

    # extract the lag 0 covariance
    R0 = Rtilde[0:D,0:D]

    kronILv = np.kron(np.eye(L),v)
    kronILw = np.kron(np.eye(L), w)
    kron1Lv = np.kron(np.ones((L,1)) , v )
    kron1Lw = np.kron(np.ones((L,1)), w)

    # covariance vector
    r = np.kron(np.eye(2*L),v).T @ np.kron(np.eye(2),blkR1) @ np.concatenate((kron1Lv,kron1Lw),axis=0)

    # covariance matrix
    _R1 = np.concatenate( ( np.kron(np.ones((1,2)),kronILv.T), np.kron(np.ones((1,2)),kronILw.T) ), axis=0 )
    _R2 = np.kron(np.eye(2), Rtilde)
    _R3 = scipy.linalg.block_diag(kronILv, kronILw)
    R = _R1 @ _R2 @ _R3

    res = (v.T @ R0 @ v) - (r.T @ np.linalg.pinv(R) @ r)

    return res

def nlcv(x):
    D=x.shape[0]
    x1 = x[0:D]
    return x1.T@x1

def nlcw(x):
    D=x.shape[0]
    x2 = x[D+1:2*D+1]
    return x2.T@x2

def gca_obj(x, blkR1, Rtilde, blkR1r, Rtilder):

    """

    NB: x = (v,w) driven then driving

    :param x: np.concatenate((v,w),axis=0), driven then driving
    :param blkR1: block covariance (the one with zero matrices, see paper)
    :param Rtilde: block covariance (the block toeplitz one, see paper)
    :param blkR1r: block covariance with time reversed
    :param Rtilder: block covariance with time reversed
    :return: log(MMSE_full) - log(MMSE_reduced)
    """

    x = np.atleast_2d(np.matrix.flatten(x)).T
    twoD = x.shape[0] # twice the dimensionality so that x = (v,w)

    #ensure that the input is even
    assert (twoD%2)==0

    # ensure that blkR1 and Rtilde have same shape
    assert blkR1.shape==Rtilde.shape
    assert blkR1r.shape==Rtilder.shape

    D = twoD//2

    # extract the spatial filters of driven and driving signals
    v = x[0:D]
    w = x[D:twoD+1]

    # deduce L, the max lag parameter
    L=blkR1.shape[0]//D

    # forward time objective
    Gf = np.log(phi_f(v,w,blkR1,Rtilde,L)) - np.log(phi_r(v,blkR1,Rtilde,L))
    #Gf = phi_f(v,w,blkR1,Rtilde,L)/phi_r(v,blkR1,Rtilde,L) - 1

    # reverse time objective
    # note that v and w are swapped here on purpose!
    Gr = np.log(phi_f(w,v,blkR1r,Rtilder,L)) - np.log(phi_r(w,blkR1r,Rtilder,L))
    #Gr = phi_f(w,v,blkR1r,Rtilder,L)/phi_r(w,blkR1r,Rtilder,L) - 1

    G = 0.5*(Gf+Gr)

    return G

def gca_obj_v(v, w, blkR1, Rtilde, blkR1r, Rtilder):
    return gca_obj(np.concatenate( (v, w) ,axis=0 ), blkR1, Rtilde, blkR1r, Rtilder)

def gca_obj_w(w, v, blkR1, Rtilde, blkR1r, Rtilder):
    return gca_obj(np.concatenate( (v, w) ,axis=0 ), blkR1, Rtilde, blkR1r, Rtilder)

def unit_norm(x):
    return x.T@x