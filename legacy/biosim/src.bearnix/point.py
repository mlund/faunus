# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _point

def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class point(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, point, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, point, name)
    def __repr__(self):
        return "<C point instance at %s>" % (self.this,)
    __swig_setmethods__["x"] = _point.point_x_set
    __swig_getmethods__["x"] = _point.point_x_get
    if _newclass:x = property(_point.point_x_get, _point.point_x_set)
    __swig_setmethods__["y"] = _point.point_y_set
    __swig_getmethods__["y"] = _point.point_y_get
    if _newclass:y = property(_point.point_y_get, _point.point_y_set)
    __swig_setmethods__["z"] = _point.point_z_set
    __swig_getmethods__["z"] = _point.point_z_get
    if _newclass:z = property(_point.point_z_get, _point.point_z_set)
    def __init__(self, *args):
        _swig_setattr(self, point, 'this', _point.new_point(*args))
        _swig_setattr(self, point, 'thisown', 1)
    def clear(*args): return _point.point_clear(*args)
    def len(*args): return _point.point_len(*args)
    def sqdist(*args): return _point.point_sqdist(*args)
    def invdist(*args): return _point.point_invdist(*args)
    def dist(*args): return _point.point_dist(*args)
    def dot(*args): return _point.point_dot(*args)
    def __neg__(*args): return _point.point___neg__(*args)
    def __mul__(*args): return _point.point___mul__(*args)
    def __add__(*args): return _point.point___add__(*args)
    def __sub__(*args): return _point.point___sub__(*args)
    def __iadd__(*args): return _point.point___iadd__(*args)
    def __del__(self, destroy=_point.delete_point):
        try:
            if self.thisown: destroy(self)
        except: pass

class pointPtr(point):
    def __init__(self, this):
        _swig_setattr(self, point, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, point, 'thisown', 0)
        _swig_setattr(self, point,self.__class__,point)
_point.point_swigregister(pointPtr)

class particle(point):
    __swig_setmethods__ = {}
    for _s in [point]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, particle, name, value)
    __swig_getmethods__ = {}
    for _s in [point]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, particle, name)
    def __repr__(self):
        return "<C particle instance at %s>" % (self.this,)
    __swig_setmethods__["charge"] = _point.particle_charge_set
    __swig_getmethods__["charge"] = _point.particle_charge_get
    if _newclass:charge = property(_point.particle_charge_get, _point.particle_charge_set)
    __swig_setmethods__["radius"] = _point.particle_radius_set
    __swig_getmethods__["radius"] = _point.particle_radius_get
    if _newclass:radius = property(_point.particle_radius_get, _point.particle_radius_set)
    __swig_setmethods__["mw"] = _point.particle_mw_set
    __swig_getmethods__["mw"] = _point.particle_mw_get
    if _newclass:mw = property(_point.particle_mw_get, _point.particle_mw_set)
    __swig_setmethods__["id"] = _point.particle_id_set
    __swig_getmethods__["id"] = _point.particle_id_get
    if _newclass:id = property(_point.particle_id_get, _point.particle_id_set)
    __swig_setmethods__["hydr"] = _point.particle_hydr_set
    __swig_getmethods__["hydr"] = _point.particle_hydr_get
    if _newclass:hydr = property(_point.particle_hydr_get, _point.particle_hydr_set)
    def __init__(self, *args):
        _swig_setattr(self, particle, 'this', _point.new_particle(*args))
        _swig_setattr(self, particle, 'thisown', 1)
    def overlap(*args): return _point.particle_overlap(*args)
    def potential(*args): return _point.particle_potential(*args)
    def __del__(self, destroy=_point.delete_particle):
        try:
            if self.thisown: destroy(self)
        except: pass

class particlePtr(particle):
    def __init__(self, this):
        _swig_setattr(self, particle, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, particle, 'thisown', 0)
        _swig_setattr(self, particle,self.__class__,particle)
_point.particle_swigregister(particlePtr)

class spherical(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, spherical, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, spherical, name)
    def __repr__(self):
        return "<C spherical instance at %s>" % (self.this,)
    __swig_setmethods__["r"] = _point.spherical_r_set
    __swig_getmethods__["r"] = _point.spherical_r_get
    if _newclass:r = property(_point.spherical_r_get, _point.spherical_r_set)
    __swig_setmethods__["theta"] = _point.spherical_theta_set
    __swig_getmethods__["theta"] = _point.spherical_theta_get
    if _newclass:theta = property(_point.spherical_theta_get, _point.spherical_theta_set)
    __swig_setmethods__["phi"] = _point.spherical_phi_set
    __swig_getmethods__["phi"] = _point.spherical_phi_get
    if _newclass:phi = property(_point.spherical_phi_get, _point.spherical_phi_set)
    def __init__(self, *args):
        _swig_setattr(self, spherical, 'this', _point.new_spherical(*args))
        _swig_setattr(self, spherical, 'thisown', 1)
    def cartesian(*args): return _point.spherical_cartesian(*args)
    def __del__(self, destroy=_point.delete_spherical):
        try:
            if self.thisown: destroy(self)
        except: pass

class sphericalPtr(spherical):
    def __init__(self, this):
        _swig_setattr(self, spherical, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, spherical, 'thisown', 0)
        _swig_setattr(self, spherical,self.__class__,spherical)
_point.spherical_swigregister(sphericalPtr)


