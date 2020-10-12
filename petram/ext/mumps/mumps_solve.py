# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.1
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _mumps_solve
else:
    import _mumps_solve

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class mumps_complex(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    r = property(_mumps_solve.mumps_complex_r_get, _mumps_solve.mumps_complex_r_set)
    i = property(_mumps_solve.mumps_complex_i_get, _mumps_solve.mumps_complex_i_set)

    def __init__(self):
        _mumps_solve.mumps_complex_swiginit(self, _mumps_solve.new_mumps_complex())
    __swig_destroy__ = _mumps_solve.delete_mumps_complex

# Register mumps_complex in _mumps_solve:
_mumps_solve.mumps_complex_swigregister(mumps_complex)

class mumps_double_complex(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    r = property(_mumps_solve.mumps_double_complex_r_get, _mumps_solve.mumps_double_complex_r_set)
    i = property(_mumps_solve.mumps_double_complex_i_get, _mumps_solve.mumps_double_complex_i_set)

    def __init__(self):
        _mumps_solve.mumps_double_complex_swiginit(self, _mumps_solve.new_mumps_double_complex())
    __swig_destroy__ = _mumps_solve.delete_mumps_double_complex

# Register mumps_double_complex in _mumps_solve:
_mumps_solve.mumps_double_complex_swigregister(mumps_double_complex)

MUMPS_ARITH_s = _mumps_solve.MUMPS_ARITH_s
MUMPS_ARITH_d = _mumps_solve.MUMPS_ARITH_d
MUMPS_ARITH_c = _mumps_solve.MUMPS_ARITH_c
MUMPS_ARITH_z = _mumps_solve.MUMPS_ARITH_z
MUMPS_ARITH_REAL = _mumps_solve.MUMPS_ARITH_REAL
MUMPS_ARITH_CMPLX = _mumps_solve.MUMPS_ARITH_CMPLX
MUMPS_ARITH_SINGLE = _mumps_solve.MUMPS_ARITH_SINGLE
MUMPS_ARITH_DBL = _mumps_solve.MUMPS_ARITH_DBL
JOB_INIT = _mumps_solve.JOB_INIT
JOB_END = _mumps_solve.JOB_END
JOB_ANALYSIS = _mumps_solve.JOB_ANALYSIS
JOB_FACTORIZATION = _mumps_solve.JOB_FACTORIZATION
JOB_SOLVE = _mumps_solve.JOB_SOLVE
JOB_1_2 = _mumps_solve.JOB_1_2
JOB_2_3 = _mumps_solve.JOB_2_3
JOB_1_2_3 = _mumps_solve.JOB_1_2_3
USE_COMM_WORLD = _mumps_solve.USE_COMM_WORLD
class DMUMPS(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        _mumps_solve.DMUMPS_swiginit(self, _mumps_solve.new_DMUMPS(*args))
    __swig_destroy__ = _mumps_solve.delete_DMUMPS

    def run(self):

        import sys
        sys.stdout.flush()
        sys.stderr.flush()  


        return _mumps_solve.DMUMPS_run(self)


    def finish(self):
        return _mumps_solve.DMUMPS_finish(self)

    def set_cntl(self, i, v):
        return _mumps_solve.DMUMPS_set_cntl(self, i, v)

    def set_icntl(self, i, v):
        return _mumps_solve.DMUMPS_set_icntl(self, i, v)

    def get_cntl(self, i):
        return _mumps_solve.DMUMPS_get_cntl(self, i)

    def get_icntl(self, i):
        return _mumps_solve.DMUMPS_get_icntl(self, i)

    def set_job(self, n):
        return _mumps_solve.DMUMPS_set_job(self, n)

    def set_n(self, n):
        return _mumps_solve.DMUMPS_set_n(self, n)

    def set_nz(self, nz):
        return _mumps_solve.DMUMPS_set_nz(self, nz)

    def set_irn(self, irn):
        return _mumps_solve.DMUMPS_set_irn(self, irn)

    def set_jcn(self, jcn):
        return _mumps_solve.DMUMPS_set_jcn(self, jcn)

    def set_nz_loc(self, nz):
        return _mumps_solve.DMUMPS_set_nz_loc(self, nz)

    def set_irn_loc(self, irn):
        return _mumps_solve.DMUMPS_set_irn_loc(self, irn)

    def set_jcn_loc(self, jcn):
        return _mumps_solve.DMUMPS_set_jcn_loc(self, jcn)

    def set_a(self, a):
        return _mumps_solve.DMUMPS_set_a(self, a)

    def set_a_loc(self, a):
        return _mumps_solve.DMUMPS_set_a_loc(self, a)

    def set_rhs(self, rhs):
        return _mumps_solve.DMUMPS_set_rhs(self, rhs)

    def set_lrhs_nrhs(self, lrhs, nrhs):
        return _mumps_solve.DMUMPS_set_lrhs_nrhs(self, lrhs, nrhs)

    def set_saveparam(self, prefix, dir):
        return _mumps_solve.DMUMPS_set_saveparam(self, prefix, dir)

    def get_rhs(self):
        return _mumps_solve.DMUMPS_get_rhs(self)

    def set_ictrl(self, i):
        return _mumps_solve.DMUMPS_set_ictrl(self, i)

    def get_struct(self):
        return _mumps_solve.DMUMPS_get_struct(self)

    def get_info(self, i):
        return _mumps_solve.DMUMPS_get_info(self, i)

    def get_infog(self, i):
        return _mumps_solve.DMUMPS_get_infog(self, i)

    def get_rinfo(self, i):
        return _mumps_solve.DMUMPS_get_rinfo(self, i)

    def get_rinfog(self, i):
        return _mumps_solve.DMUMPS_get_rinfog(self, i)

    def get_version_number(self):
        return _mumps_solve.DMUMPS_get_version_number(self)

    def get_real_rhs(self):
        return _mumps_solve.DMUMPS_get_real_rhs(self)

# Register DMUMPS in _mumps_solve:
_mumps_solve.DMUMPS_swigregister(DMUMPS)

class ZMUMPS(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        _mumps_solve.ZMUMPS_swiginit(self, _mumps_solve.new_ZMUMPS(*args))
    __swig_destroy__ = _mumps_solve.delete_ZMUMPS

    def run(self):

        import sys
        sys.stdout.flush()
        sys.stderr.flush()  


        return _mumps_solve.ZMUMPS_run(self)


    def finish(self):
        return _mumps_solve.ZMUMPS_finish(self)

    def set_cntl(self, i, v):
        return _mumps_solve.ZMUMPS_set_cntl(self, i, v)

    def set_icntl(self, i, v):
        return _mumps_solve.ZMUMPS_set_icntl(self, i, v)

    def get_cntl(self, i):
        return _mumps_solve.ZMUMPS_get_cntl(self, i)

    def get_icntl(self, i):
        return _mumps_solve.ZMUMPS_get_icntl(self, i)

    def set_job(self, n):
        return _mumps_solve.ZMUMPS_set_job(self, n)

    def set_n(self, n):
        return _mumps_solve.ZMUMPS_set_n(self, n)

    def set_nz(self, nz):
        return _mumps_solve.ZMUMPS_set_nz(self, nz)

    def set_irn(self, irn):
        return _mumps_solve.ZMUMPS_set_irn(self, irn)

    def set_jcn(self, jcn):
        return _mumps_solve.ZMUMPS_set_jcn(self, jcn)

    def set_nz_loc(self, nz):
        return _mumps_solve.ZMUMPS_set_nz_loc(self, nz)

    def set_irn_loc(self, irn):
        return _mumps_solve.ZMUMPS_set_irn_loc(self, irn)

    def set_jcn_loc(self, jcn):
        return _mumps_solve.ZMUMPS_set_jcn_loc(self, jcn)

    def set_a(self, a):
        return _mumps_solve.ZMUMPS_set_a(self, a)

    def set_a_loc(self, a):
        return _mumps_solve.ZMUMPS_set_a_loc(self, a)

    def set_rhs(self, rhs):
        return _mumps_solve.ZMUMPS_set_rhs(self, rhs)

    def set_lrhs_nrhs(self, lrhs, nrhs):
        return _mumps_solve.ZMUMPS_set_lrhs_nrhs(self, lrhs, nrhs)

    def set_saveparam(self, prefix, dir):
        return _mumps_solve.ZMUMPS_set_saveparam(self, prefix, dir)

    def get_rhs(self):
        return _mumps_solve.ZMUMPS_get_rhs(self)

    def set_ictrl(self, i):
        return _mumps_solve.ZMUMPS_set_ictrl(self, i)

    def get_struct(self):
        return _mumps_solve.ZMUMPS_get_struct(self)

    def get_info(self, i):
        return _mumps_solve.ZMUMPS_get_info(self, i)

    def get_infog(self, i):
        return _mumps_solve.ZMUMPS_get_infog(self, i)

    def get_rinfo(self, i):
        return _mumps_solve.ZMUMPS_get_rinfo(self, i)

    def get_rinfog(self, i):
        return _mumps_solve.ZMUMPS_get_rinfog(self, i)

    def get_version_number(self):
        return _mumps_solve.ZMUMPS_get_version_number(self)

    def get_real_rhs(self):
        return _mumps_solve.ZMUMPS_get_real_rhs(self)

    def get_imag_rhs(self):
        return _mumps_solve.ZMUMPS_get_imag_rhs(self)

# Register ZMUMPS in _mumps_solve:
_mumps_solve.ZMUMPS_swigregister(ZMUMPS)

class SMUMPS(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        _mumps_solve.SMUMPS_swiginit(self, _mumps_solve.new_SMUMPS(*args))
    __swig_destroy__ = _mumps_solve.delete_SMUMPS

    def run(self):

        import sys
        sys.stdout.flush()
        sys.stderr.flush()  


        return _mumps_solve.SMUMPS_run(self)


    def finish(self):
        return _mumps_solve.SMUMPS_finish(self)

    def set_cntl(self, i, v):
        return _mumps_solve.SMUMPS_set_cntl(self, i, v)

    def set_icntl(self, i, v):
        return _mumps_solve.SMUMPS_set_icntl(self, i, v)

    def get_cntl(self, i):
        return _mumps_solve.SMUMPS_get_cntl(self, i)

    def get_icntl(self, i):
        return _mumps_solve.SMUMPS_get_icntl(self, i)

    def set_job(self, n):
        return _mumps_solve.SMUMPS_set_job(self, n)

    def set_n(self, n):
        return _mumps_solve.SMUMPS_set_n(self, n)

    def set_nz(self, nz):
        return _mumps_solve.SMUMPS_set_nz(self, nz)

    def set_irn(self, irn):
        return _mumps_solve.SMUMPS_set_irn(self, irn)

    def set_jcn(self, jcn):
        return _mumps_solve.SMUMPS_set_jcn(self, jcn)

    def set_nz_loc(self, nz):
        return _mumps_solve.SMUMPS_set_nz_loc(self, nz)

    def set_irn_loc(self, irn):
        return _mumps_solve.SMUMPS_set_irn_loc(self, irn)

    def set_jcn_loc(self, jcn):
        return _mumps_solve.SMUMPS_set_jcn_loc(self, jcn)

    def set_a(self, a):
        return _mumps_solve.SMUMPS_set_a(self, a)

    def set_a_loc(self, a):
        return _mumps_solve.SMUMPS_set_a_loc(self, a)

    def set_rhs(self, rhs):
        return _mumps_solve.SMUMPS_set_rhs(self, rhs)

    def set_lrhs_nrhs(self, lrhs, nrhs):
        return _mumps_solve.SMUMPS_set_lrhs_nrhs(self, lrhs, nrhs)

    def set_saveparam(self, prefix, dir):
        return _mumps_solve.SMUMPS_set_saveparam(self, prefix, dir)

    def get_rhs(self):
        return _mumps_solve.SMUMPS_get_rhs(self)

    def set_ictrl(self, i):
        return _mumps_solve.SMUMPS_set_ictrl(self, i)

    def get_struct(self):
        return _mumps_solve.SMUMPS_get_struct(self)

    def get_info(self, i):
        return _mumps_solve.SMUMPS_get_info(self, i)

    def get_infog(self, i):
        return _mumps_solve.SMUMPS_get_infog(self, i)

    def get_rinfo(self, i):
        return _mumps_solve.SMUMPS_get_rinfo(self, i)

    def get_rinfog(self, i):
        return _mumps_solve.SMUMPS_get_rinfog(self, i)

    def get_version_number(self):
        return _mumps_solve.SMUMPS_get_version_number(self)

    def get_real_rhs(self):
        return _mumps_solve.SMUMPS_get_real_rhs(self)

# Register SMUMPS in _mumps_solve:
_mumps_solve.SMUMPS_swigregister(SMUMPS)

class CMUMPS(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        _mumps_solve.CMUMPS_swiginit(self, _mumps_solve.new_CMUMPS(*args))
    __swig_destroy__ = _mumps_solve.delete_CMUMPS

    def run(self):

        import sys
        sys.stdout.flush()
        sys.stderr.flush()  


        return _mumps_solve.CMUMPS_run(self)


    def finish(self):
        return _mumps_solve.CMUMPS_finish(self)

    def set_cntl(self, i, v):
        return _mumps_solve.CMUMPS_set_cntl(self, i, v)

    def set_icntl(self, i, v):
        return _mumps_solve.CMUMPS_set_icntl(self, i, v)

    def get_cntl(self, i):
        return _mumps_solve.CMUMPS_get_cntl(self, i)

    def get_icntl(self, i):
        return _mumps_solve.CMUMPS_get_icntl(self, i)

    def set_job(self, n):
        return _mumps_solve.CMUMPS_set_job(self, n)

    def set_n(self, n):
        return _mumps_solve.CMUMPS_set_n(self, n)

    def set_nz(self, nz):
        return _mumps_solve.CMUMPS_set_nz(self, nz)

    def set_irn(self, irn):
        return _mumps_solve.CMUMPS_set_irn(self, irn)

    def set_jcn(self, jcn):
        return _mumps_solve.CMUMPS_set_jcn(self, jcn)

    def set_nz_loc(self, nz):
        return _mumps_solve.CMUMPS_set_nz_loc(self, nz)

    def set_irn_loc(self, irn):
        return _mumps_solve.CMUMPS_set_irn_loc(self, irn)

    def set_jcn_loc(self, jcn):
        return _mumps_solve.CMUMPS_set_jcn_loc(self, jcn)

    def set_a(self, a):
        return _mumps_solve.CMUMPS_set_a(self, a)

    def set_a_loc(self, a):
        return _mumps_solve.CMUMPS_set_a_loc(self, a)

    def set_rhs(self, rhs):
        return _mumps_solve.CMUMPS_set_rhs(self, rhs)

    def set_lrhs_nrhs(self, lrhs, nrhs):
        return _mumps_solve.CMUMPS_set_lrhs_nrhs(self, lrhs, nrhs)

    def set_saveparam(self, prefix, dir):
        return _mumps_solve.CMUMPS_set_saveparam(self, prefix, dir)

    def get_rhs(self):
        return _mumps_solve.CMUMPS_get_rhs(self)

    def set_ictrl(self, i):
        return _mumps_solve.CMUMPS_set_ictrl(self, i)

    def get_struct(self):
        return _mumps_solve.CMUMPS_get_struct(self)

    def get_info(self, i):
        return _mumps_solve.CMUMPS_get_info(self, i)

    def get_infog(self, i):
        return _mumps_solve.CMUMPS_get_infog(self, i)

    def get_rinfo(self, i):
        return _mumps_solve.CMUMPS_get_rinfo(self, i)

    def get_rinfog(self, i):
        return _mumps_solve.CMUMPS_get_rinfog(self, i)

    def get_version_number(self):
        return _mumps_solve.CMUMPS_get_version_number(self)

    def get_real_rhs(self):
        return _mumps_solve.CMUMPS_get_real_rhs(self)

    def get_imag_rhs(self):
        return _mumps_solve.CMUMPS_get_imag_rhs(self)

# Register CMUMPS in _mumps_solve:
_mumps_solve.CMUMPS_swigregister(CMUMPS)


def libmumps_solve_example_dist(comm):
    return _mumps_solve.libmumps_solve_example_dist(comm)

def libmumps_solve_example(comm):
    return _mumps_solve.libmumps_solve_example(comm)

def SIZEOF_MUMPS_INT():
    return _mumps_solve.SIZEOF_MUMPS_INT()

def i_array(arr):
    return _mumps_solve.i_array(arr)

def s_array(arr):
    return _mumps_solve.s_array(arr)

def d_array(arr):
    return _mumps_solve.d_array(arr)

def c_array(arr):
    return _mumps_solve.c_array(arr)

def z_array(arr):
    return _mumps_solve.z_array(arr)

def c_real_array(arr):
    return _mumps_solve.c_real_array(arr)

def z_real_array(arr):
    return _mumps_solve.z_real_array(arr)

def i_array_getitem(arr, i):
    return _mumps_solve.i_array_getitem(arr, i)

def d_array_getitem(arr, i):
    return _mumps_solve.d_array_getitem(arr, i)

def s_array_getitem(arr, i):
    return _mumps_solve.s_array_getitem(arr, i)

def z_array_real_getitem(arr, i):
    return _mumps_solve.z_array_real_getitem(arr, i)

def c_array_real_getitem(arr, i):
    return _mumps_solve.c_array_real_getitem(arr, i)

def z_array_imag_getitem(arr, i):
    return _mumps_solve.z_array_imag_getitem(arr, i)

def c_array_imag_getitem(arr, i):
    return _mumps_solve.c_array_imag_getitem(arr, i)

def i_to_list(A, l):
    return [i_array_getitem(A, i) for i in range(l)]
def d_to_list(A, l):
    return [d_array_getitem(A, i) for i in range(l)]
def s_to_list(A, l):
    return [s_array_getitem(A, i) for i in range(l)]
def c_to_list(A, l):
    return [c_array_real_getitem(A, i) +
	    1j* c_array_imag_getitem(A, i) for i in range(l)]
def z_to_list(A, l):
    return [z_array_real_getitem(A, i) +
	    1j* z_array_imag_getitem(A, i) for i in range(l)]        



