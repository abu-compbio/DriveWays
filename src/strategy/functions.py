# This is
import os
import numpy as np
import cProfile, pstats, io
def cov_over_mutex(module, t=None):

    return np.float128(module.cov)/ np.float128(module.mutex) if t== None else (np.float128(module.cov)/ np.float128(module.mutex))*np.float128(t)

def cov(module, t=None):
    return np.float128(module.cov) if t== None else np.float128(module.cov)/np.float128(t)

def cov_mutx(module, t=None):
    return np.float128(module.cov)* np.float128(module.mutex) if t== None else (np.float128(module.cov)* np.float128(module.mutex))*np.float128(t)

def cov_cov_mutex(module, t=None):
    # This is similar to the buggy implementation in java, discovered on 08/02/2019
    # Not needed most probably
    return np.float128(module.cov)/ np.float128(module.mutex)  if t== None else np.float128(module.cov)/np.float128(t)


fncs = {
        'cov/mutex': cov_over_mutex,
        'cov': cov,
        'cov*mutex': cov_mutx,
        'cov>cov/mutex': cov_cov_mutex
}


def mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)





# This a decorator that help in profiling the code, and get
def profile(fnc):

    """A decorator that uses cProfile to profile a function"""

    def inner(*args, **kwargs):

        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner
