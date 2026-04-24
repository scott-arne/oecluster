/* -*- mode: c++ -*- */
// swig/oecluster.i
// SWIG interface file for oecluster Python bindings
%module _oecluster

%{
#include "oecluster/oecluster.h"
#include "oecluster/Error.h"
#include "oecluster/DistanceMetric.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/ThreadPool.h"
#include "oecluster/PDist.h"
#include "oecluster/CDist.h"
#include "oecluster/DistanceMatrix.h"
#include "oecluster/metrics/FingerprintMetric.h"
#include "oecluster/metrics/ROCSMetric.h"
#include "oecluster/metrics/SuperposeMetric.h"

#include <oechem.h>
#include <oebio.h>
#include <oespruce.h>
#include <oegrid.h>

using namespace OECluster;
%}

// ============================================================================
// Forward declarations for cross-module SWIG type resolution
// ============================================================================
// These enable typemaps for OpenEye types whose definitions live in the
// OpenEye SWIG runtime (v4). Only types you actually use in your wrapped API
// need full #include — forward declarations suffice for the typemaps.

namespace OEChem {
    class OEMolBase;
    class OEMCMolBase;
    class OEMol;
    class OEGraphMol;
    class OEAtomBase;
    class OEBondBase;
    class OEConfBase;
    class OEMatchBase;
    class OEMolDatabase;
    class oemolistream;
    class oemolostream;
    class OEQMol;
    class OEResidue;
    class OEUniMolecularRxn;
}

namespace OEBio {
    class OEDesignUnit;
    class OEHierView;
    class OEHierResidue;
    class OEHierFragment;
    class OEHierChain;
    class OEInteractionHint;
    class OEInteractionHintContainer;
}

namespace OEDocking {
    class OEReceptor;
}

namespace OEPlatform {
    class oeifstream;
    class oeofstream;
    class oeisstream;
    class oeosstream;
}

namespace OESystem {
    class OEScalarGrid;
    class OERecord;
    class OEMolRecord;
}

// ============================================================================
// Cross-runtime SWIG compatibility layer
// ============================================================================
// OpenEye's Python bindings use SWIG runtime v4; our module uses v5.
// Since the runtimes are separate, SWIG_TypeQuery cannot access OpenEye types.
// We use Python isinstance for type safety and directly extract the void*
// pointer from the SwigPyObject struct layout (stable across SWIG versions).
//
// This approach enables passing OpenEye objects between Python and C++ without
// serialization. The macros below generate the boilerplate for each type.

%{
// Minimal SwigPyObject layout compatible across SWIG runtime versions.
// The actual struct may have more fields, but ptr is always first after
// PyObject_HEAD.
struct _SwigPyObjectCompat {
    PyObject_HEAD
    void *ptr;
};

static void* _oecluster_extract_swig_ptr(PyObject* obj) {
    PyObject* thisAttr = PyObject_GetAttrString(obj, "this");
    if (!thisAttr) {
        PyErr_Clear();
        return NULL;
    }
    void* ptr = ((_SwigPyObjectCompat*)thisAttr)->ptr;
    Py_DECREF(thisAttr);
    return ptr;
}

// ---- Type checker generator macro ----
// Generates a cached isinstance checker for an OpenEye Python type.
// TAG:    identifier suffix (e.g., oemolbase)
// MODULE: Python module string (e.g., "openeye.oechem")
// CLASS:  Python class name string (e.g., "OEMolBase")
#define DEFINE_OE_TYPE_CHECKER(TAG, MODULE, CLASS) \
    static PyObject* _oecluster_oe_##TAG##_type = NULL; \
    static bool _oecluster_is_##TAG(PyObject* obj) { \
        if (!_oecluster_oe_##TAG##_type) { \
            PyObject* mod = PyImport_ImportModule(MODULE); \
            if (mod) { \
                _oecluster_oe_##TAG##_type = PyObject_GetAttrString(mod, CLASS); \
                Py_DECREF(mod); \
            } \
            if (!_oecluster_oe_##TAG##_type) return false; \
        } \
        return PyObject_IsInstance(obj, _oecluster_oe_##TAG##_type) == 1; \
    }

// ---- Molecule types (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oemolbase,    "openeye.oechem", "OEMolBase")
DEFINE_OE_TYPE_CHECKER(oemcmolbase,  "openeye.oechem", "OEMCMolBase")
DEFINE_OE_TYPE_CHECKER(oemol,        "openeye.oechem", "OEMol")
DEFINE_OE_TYPE_CHECKER(oegraphmol,   "openeye.oechem", "OEGraphMol")
DEFINE_OE_TYPE_CHECKER(oeqmol,       "openeye.oechem", "OEQMol")

// ---- Atom / bond / conformer / residue (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeatombase,   "openeye.oechem", "OEAtomBase")
DEFINE_OE_TYPE_CHECKER(oebondbase,   "openeye.oechem", "OEBondBase")
DEFINE_OE_TYPE_CHECKER(oeconfbase,   "openeye.oechem", "OEConfBase")
DEFINE_OE_TYPE_CHECKER(oeresidue,    "openeye.oechem", "OEResidue")
DEFINE_OE_TYPE_CHECKER(oematchbase,  "openeye.oechem", "OEMatchBase")

// ---- Molecule I/O (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oemolistream, "openeye.oechem", "oemolistream")
DEFINE_OE_TYPE_CHECKER(oemolostream, "openeye.oechem", "oemolostream")
DEFINE_OE_TYPE_CHECKER(oemoldatabase,"openeye.oechem", "OEMolDatabase")

// ---- Reactions (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeunimolecularrxn, "openeye.oechem", "OEUniMolecularRxn")

// ---- Platform streams (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeifstream,   "openeye.oechem", "oeifstream")
DEFINE_OE_TYPE_CHECKER(oeofstream,   "openeye.oechem", "oeofstream")
DEFINE_OE_TYPE_CHECKER(oeisstream,   "openeye.oechem", "oeisstream")
DEFINE_OE_TYPE_CHECKER(oeosstream,   "openeye.oechem", "oeosstream")

// ---- Records (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oerecord,     "openeye.oechem", "OERecord")
DEFINE_OE_TYPE_CHECKER(oemolrecord,  "openeye.oechem", "OEMolRecord")

// ---- Bio / hierarchy (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oedesignunit, "openeye.oechem", "OEDesignUnit")
DEFINE_OE_TYPE_CHECKER(oehierview,   "openeye.oechem", "OEHierView")
DEFINE_OE_TYPE_CHECKER(oehierresidue,"openeye.oechem", "OEHierResidue")
DEFINE_OE_TYPE_CHECKER(oehierfragment,"openeye.oechem","OEHierFragment")
DEFINE_OE_TYPE_CHECKER(oehierchain,  "openeye.oechem", "OEHierChain")
DEFINE_OE_TYPE_CHECKER(oeinteractionhint,          "openeye.oechem", "OEInteractionHint")
DEFINE_OE_TYPE_CHECKER(oeinteractionhintcontainer, "openeye.oechem", "OEInteractionHintContainer")

// ---- Grid (openeye.oegrid) ----
DEFINE_OE_TYPE_CHECKER(oescalargrid, "openeye.oegrid", "OEScalarGrid")

// ---- Docking (openeye.oedocking) ----
DEFINE_OE_TYPE_CHECKER(oereceptor,   "openeye.oedocking", "OEReceptor")

#undef DEFINE_OE_TYPE_CHECKER

// ---- OEScalarGrid return-type helper (zero-copy pointer swap) ----
static PyObject* _oecluster_wrap_as_oe_grid(OESystem::OEScalarGrid* grid) {
    if (!grid) {
        Py_RETURN_NONE;
    }
    PyObject* oegrid_mod = PyImport_ImportModule("openeye.oegrid");
    if (!oegrid_mod) {
        delete grid;
        return NULL;
    }
    PyObject* grid_cls = PyObject_GetAttrString(oegrid_mod, "OEScalarGrid");
    Py_DECREF(oegrid_mod);
    if (!grid_cls) {
        delete grid;
        return NULL;
    }
    PyObject* oe_grid = PyObject_CallNoArgs(grid_cls);
    Py_DECREF(grid_cls);
    if (!oe_grid) {
        delete grid;
        return NULL;
    }
    PyObject* thisAttr = PyObject_GetAttrString(oe_grid, "this");
    if (!thisAttr) {
        PyErr_Clear();
        Py_DECREF(oe_grid);
        delete grid;
        return NULL;
    }
    _SwigPyObjectCompat* swig_this = (_SwigPyObjectCompat*)thisAttr;
    delete reinterpret_cast<OESystem::OEScalarGrid*>(swig_this->ptr);
    swig_this->ptr = grid;
    Py_DECREF(thisAttr);
    return oe_grid;
}
%}

// ============================================================================
// Typemap generator macros
// ============================================================================

// Generate const-ref and non-const-ref typemaps for a cross-runtime OpenEye type.
// CPP_TYPE: fully qualified C++ type (e.g., OEChem::OEMolBase)
// CHECKER:  isinstance checker function name
// ERR_MSG:  error message on type mismatch
%define OE_CROSS_RUNTIME_REF_TYPEMAPS(CPP_TYPE, CHECKER, ERR_MSG)

%typemap(in) const CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (CHECKER($input)) {
            argp = _oecluster_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) const CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
}

%typemap(in) CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (CHECKER($input)) {
            argp = _oecluster_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
}

%enddef

// Generate nullable-pointer typemaps (accepts None) for a cross-runtime type.
%define OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(CPP_TYPE, CHECKER, ERR_MSG)

%typemap(in) const CPP_TYPE* (void *argp = 0, int res = 0) {
    if ($input == Py_None) {
        $1 = NULL;
    } else {
        res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
        if (!SWIG_IsOK(res)) {
            if (CHECKER($input)) {
                argp = _oecluster_extract_swig_ptr($input);
                if (argp) res = SWIG_OK;
            }
        }
        if (!SWIG_IsOK(res)) {
            SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
        }
        $1 = reinterpret_cast< $1_ltype >(argp);
    }
}

%typemap(typecheck, precedence=10) const CPP_TYPE* {
    if ($input == Py_None) {
        $1 = 1;
    } else {
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $descriptor, 0);
        $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
    }
}

%enddef

// ============================================================================
// Typemap declarations for all OpenEye types
// ============================================================================
// Each type gets const-ref and non-const-ref typemaps. Types that commonly
// appear as optional parameters also get nullable-pointer typemaps.
// These are inert until a wrapped function signature uses the type.

// ---- Molecule hierarchy (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolBase,    _oecluster_is_oemolbase,    "Expected OEMolBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMCMolBase,  _oecluster_is_oemcmolbase,  "Expected OEMCMolBase-derived object (OEMCMolBase or OEMol).")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMol,        _oecluster_is_oemol,        "Expected OEMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEGraphMol,   _oecluster_is_oegraphmol,   "Expected OEGraphMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEQMol,       _oecluster_is_oeqmol,       "Expected OEQMol object.")

// ---- Atom / bond / conformer / residue / match (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEAtomBase,   _oecluster_is_oeatombase,   "Expected OEAtomBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEBondBase,   _oecluster_is_oebondbase,   "Expected OEBondBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEConfBase,   _oecluster_is_oeconfbase,   "Expected OEConfBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEResidue,    _oecluster_is_oeresidue,    "Expected OEResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMatchBase,  _oecluster_is_oematchbase,  "Expected OEMatchBase-derived object.")

// ---- Molecule I/O (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolistream,  _oecluster_is_oemolistream, "Expected oemolistream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolostream,  _oecluster_is_oemolostream, "Expected oemolostream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolDatabase, _oecluster_is_oemoldatabase,"Expected OEMolDatabase object.")

// ---- Reactions (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEUniMolecularRxn, _oecluster_is_oeunimolecularrxn, "Expected OEUniMolecularRxn object.")

// ---- Platform streams (OEPlatform) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeifstream, _oecluster_is_oeifstream, "Expected oeifstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeofstream, _oecluster_is_oeofstream, "Expected oeofstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeisstream, _oecluster_is_oeisstream, "Expected oeisstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeosstream, _oecluster_is_oeosstream, "Expected oeosstream object.")

// ---- Records (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OERecord,    _oecluster_is_oerecord,    "Expected OERecord object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OEMolRecord, _oecluster_is_oemolrecord, "Expected OEMolRecord object.")

// ---- Bio / hierarchy (OEBio) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEDesignUnit,   _oecluster_is_oedesignunit, "Expected OEDesignUnit object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierView,     _oecluster_is_oehierview,   "Expected OEHierView object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierResidue,  _oecluster_is_oehierresidue,"Expected OEHierResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierFragment,  _oecluster_is_oehierfragment,"Expected OEHierFragment object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierChain,    _oecluster_is_oehierchain,  "Expected OEHierChain object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHint,          _oecluster_is_oeinteractionhint,          "Expected OEInteractionHint object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHintContainer, _oecluster_is_oeinteractionhintcontainer, "Expected OEInteractionHintContainer object.")

// ---- Grid (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OEScalarGrid, _oecluster_is_oescalargrid, "Expected OEScalarGrid-derived object.")
OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(OESystem::OEScalarGrid, _oecluster_is_oescalargrid, "Expected OEScalarGrid or None.")

// OEScalarGrid return-type typemap (wraps C++ grid as native openeye.oegrid object)
%typemap(out) OESystem::OEScalarGrid* {
    $result = _oecluster_wrap_as_oe_grid($1);
    if (!$result) SWIG_fail;
}

// ---- Docking (OEDocking) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEDocking::OEReceptor, _oecluster_is_oereceptor, "Expected OEReceptor object.")

// ============================================================================
// Domain-specific typemaps for oecluster's container-based APIs
// ============================================================================
// FingerprintMetric and SuperposeMetric (mol overload) accept lists of
// OEMolBase pointers. ROCSMetric and SuperposeMetric (DU overload) accept
// lists of shared_ptr to OEMol / OEDesignUnit — these require a copy because
// we do not own the Python-owned objects.

// ---- const std::vector<OEChem::OEMolBase*>& ----
%typemap(in) const std::vector<OEChem::OEMolBase*>& (std::vector<OEChem::OEMolBase*> temp) {
    if (!PyList_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError, "Expected a list of OEMolBase objects");
    }
    Py_ssize_t size = PyList_Size($input);
    temp.reserve(size);
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyList_GetItem($input, i);
        if (!_oecluster_is_oemolbase(item)) {
            SWIG_exception_fail(SWIG_TypeError, "List item is not an OEMolBase object");
        }
        void* ptr = _oecluster_extract_swig_ptr(item);
        if (!ptr) {
            SWIG_exception_fail(SWIG_TypeError, "List item is not an OEMolBase object");
        }
        temp.push_back(reinterpret_cast<OEChem::OEMolBase*>(ptr));
    }
    $1 = &temp;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) const std::vector<OEChem::OEMolBase*>& {
    if (PyList_Check($input) && PyList_Size($input) > 0) {
        $1 = _oecluster_is_oemolbase(PyList_GetItem($input, 0)) ? 1 : 0;
    } else {
        $1 = PyList_Check($input) ? 1 : 0;
    }
}

// ---- const std::vector<std::shared_ptr<OEChem::OEMol>>& ----
%typemap(in) const std::vector<std::shared_ptr<OEChem::OEMol>>& (std::vector<std::shared_ptr<OEChem::OEMol>> temp) {
    if (!PyList_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError, "Expected a list of OEMol objects");
    }
    Py_ssize_t size = PyList_Size($input);
    temp.reserve(size);
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyList_GetItem($input, i);
        if (!_oecluster_is_oemol(item)) {
            SWIG_exception_fail(SWIG_TypeError, "List item is not an OEMol object");
        }
        void* ptr = _oecluster_extract_swig_ptr(item);
        if (!ptr) {
            SWIG_exception_fail(SWIG_TypeError, "List item is not an OEMol object");
        }
        // Make a COPY wrapped in shared_ptr (we do not own the Python object)
        auto* mol = reinterpret_cast<OEChem::OEMol*>(ptr);
        temp.push_back(std::make_shared<OEChem::OEMol>(*mol));
    }
    $1 = &temp;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) const std::vector<std::shared_ptr<OEChem::OEMol>>& {
    if (PyList_Check($input) && PyList_Size($input) > 0) {
        $1 = _oecluster_is_oemol(PyList_GetItem($input, 0)) ? 1 : 0;
    } else {
        $1 = PyList_Check($input) ? 1 : 0;
    }
}

// ---- const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& ----
%typemap(in) const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& (std::vector<std::shared_ptr<OEBio::OEDesignUnit>> temp) {
    if (!PyList_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError, "Expected a list of OEDesignUnit objects");
    }
    Py_ssize_t size = PyList_Size($input);
    temp.reserve(size);
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyList_GetItem($input, i);
        if (!_oecluster_is_oedesignunit(item)) {
            SWIG_exception_fail(SWIG_TypeError, "List item is not an OEDesignUnit object");
        }
        void* ptr = _oecluster_extract_swig_ptr(item);
        if (!ptr) {
            SWIG_exception_fail(SWIG_TypeError, "List item is not an OEDesignUnit object");
        }
        // Make a COPY wrapped in shared_ptr (we do not own the Python object)
        auto* du = reinterpret_cast<OEBio::OEDesignUnit*>(ptr);
        temp.push_back(std::make_shared<OEBio::OEDesignUnit>(*du));
    }
    $1 = &temp;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& {
    if (PyList_Check($input) && PyList_Size($input) > 0) {
        $1 = _oecluster_is_oedesignunit(PyList_GetItem($input, 0)) ? 1 : 0;
    } else {
        $1 = PyList_Check($input) ? 1 : 0;
    }
}

// ---- std::function<void(size_t, size_t)> progress callback ----
%typemap(in) std::function<void(size_t, size_t)> {
    if ($input == Py_None) {
        $1 = nullptr;
    } else if (PyCallable_Check($input)) {
        PyObject* callback = $input;
        Py_INCREF(callback);
        // Use shared_ptr guard to ensure Py_DECREF when lambda is destroyed
        auto guard = std::shared_ptr<PyObject>(callback, [](PyObject* p) {
            PyGILState_STATE gstate = PyGILState_Ensure();
            Py_DECREF(p);
            PyGILState_Release(gstate);
        });
        $1 = [guard](size_t completed, size_t total) {
            PyGILState_STATE gstate = PyGILState_Ensure();
            PyObject* result = PyObject_CallFunction(guard.get(), "nn",
                (Py_ssize_t)completed, (Py_ssize_t)total);
            Py_XDECREF(result);
            if (PyErr_Occurred()) PyErr_Clear();
            PyGILState_Release(gstate);
        };
    } else {
        SWIG_exception_fail(SWIG_TypeError, "Expected callable or None for progress callback");
    }
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) std::function<void(size_t, size_t)> {
    $1 = ($input == Py_None || PyCallable_Check($input)) ? 1 : 0;
}

%typemap(freearg) std::function<void(size_t, size_t)> {
    // Prevent copy-on-scope-exit; the shared_ptr guard handles Py_DECREF
}

// ============================================================================
// Include STL typemaps
// ============================================================================
%include "std_string.i"
%include "std_vector.i"
%include "stdint.i"
%include "exception.i"

// Instantiate vector templates used by the Python layer
%template(StringVector) std::vector<std::string>;

// ============================================================================
// Exception handling -- general
// ============================================================================
%exception {
    try {
        $action
    } catch (const OECluster::OEClusterError& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception");
    }
}

// ============================================================================
// GIL release for pdist/cdist (override the general handler above)
// ============================================================================
%exception OECluster::pdist {
    Py_BEGIN_ALLOW_THREADS
    try {
        $action
    } catch (const OECluster::OEClusterError& e) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (const std::exception& e) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception in pdist");
    }
    Py_END_ALLOW_THREADS
}

%exception OECluster::cdist {
    Py_BEGIN_ALLOW_THREADS
    try {
        $action
    } catch (const OECluster::OEClusterError& e) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (const std::exception& e) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception in cdist");
    }
    Py_END_ALLOW_THREADS
}

// ============================================================================
// Ignore problematic members before %include
//
// StorageBackend.h pulls in <shared_mutex>, <unordered_map>, <thread>
// which SWIG cannot parse. We redeclare the classes manually instead.
//
// Also ignore private clone constructors for metrics.
// ============================================================================

// Suppress default constructors for types that don't have them
%nodefaultctor OECluster::OEClusterError;
%nodefaultctor OECluster::MetricError;
%nodefaultctor OECluster::StorageError;
%nodefaultctor OECluster::DistanceMatrix;

// Ignore Data() -- we expose _data_ptr() via %extend instead
%ignore OECluster::StorageBackend::Data;
%ignore OECluster::DenseStorage::Data;
%ignore OECluster::MMapStorage::Data;
%ignore OECluster::SparseStorage::Data;

// Ignore copy/move constructors and operators that SWIG cannot handle
%ignore OECluster::MMapStorage::MMapStorage(const MMapStorage&);
%ignore OECluster::MMapStorage::operator=;

// Ignore DistanceMetric::Clone -- returns unique_ptr, expose via %extend
%ignore OECluster::DistanceMetric::Clone;

// Ignore private constructors for metric clones
%ignore OECluster::FingerprintMetric::FingerprintMetric(std::shared_ptr<const Impl>);

// Ignore static mask parsing methods (not needed from Python)
%ignore OECluster::FingerprintMetric::ParseAtomTypeMask;
%ignore OECluster::FingerprintMetric::ParseBondTypeMask;
%ignore OECluster::ROCSMetric::ROCSMetric(std::shared_ptr<const SharedData>, const Options&);
%ignore OECluster::SuperposeMetric::SuperposeMetric(std::shared_ptr<const SharedData>, const Options&);
// Ignore SparseStorage internals that use unordered_map/shared_mutex/thread
%ignore OECluster::SparseStorage::Entries;

// ============================================================================
// Redeclare StorageBackend hierarchy for SWIG
// (instead of %include "oecluster/StorageBackend.h" which has complex STL)
// ============================================================================
namespace OECluster {

class StorageBackend {
public:
    virtual ~StorageBackend() = default;
    virtual void Set(size_t i, size_t j, double value) = 0;
    virtual double Get(size_t i, size_t j) const = 0;
    virtual size_t NumItems() const = 0;
    virtual size_t NumPairs() const = 0;
    virtual void Finalize();
};

class DenseStorage : public StorageBackend {
public:
    explicit DenseStorage(size_t n);
    void Set(size_t i, size_t j, double value) override;
    double Get(size_t i, size_t j) const override;
    size_t NumItems() const override;
    size_t NumPairs() const override;
};

class MMapStorage : public StorageBackend {
public:
    MMapStorage(const std::string& path, size_t n);
    ~MMapStorage() override;
    void Set(size_t i, size_t j, double value) override;
    double Get(size_t i, size_t j) const override;
    size_t NumItems() const override;
    size_t NumPairs() const override;
};

class SparseStorage : public StorageBackend {
public:
    SparseStorage(size_t n, double cutoff);
    ~SparseStorage() override;
    void Set(size_t i, size_t j, double value) override;
    double Get(size_t i, size_t j) const override;
    size_t NumItems() const override;
    size_t NumPairs() const override;
    void Finalize() override;
};

}  // namespace OECluster

// ============================================================================
// %extend blocks for StorageBackend and SparseStorage
// ============================================================================
%extend OECluster::StorageBackend {
    /** Return raw data pointer as integer address for numpy interop. */
    size_t _data_ptr() const {
        return reinterpret_cast<size_t>($self->Data());
    }
}

%extend OECluster::SparseStorage {
    /** Return entries as a Python list of (i, j, value) tuples. */
    PyObject* _entries() const {
        const auto& entries = $self->Entries();
        PyObject* list = PyList_New(entries.size());
        for (size_t i = 0; i < entries.size(); ++i) {
            const auto& [row, col, val] = entries[i];
            PyList_SetItem(list, i, Py_BuildValue("(nnf)",
                (Py_ssize_t)row, (Py_ssize_t)col, val));
        }
        return list;
    }
}

// ============================================================================
// Error hierarchy
// ============================================================================
%include "oecluster/Error.h"

// ============================================================================
// DistanceMetric base class
// ============================================================================
%include "oecluster/DistanceMetric.h"

// ============================================================================
// ThreadPool (PIMPL -- only expose public interface)
// Redeclare to avoid pulling in complex threading headers.
// ============================================================================
namespace OECluster {

class ThreadPool {
public:
    explicit ThreadPool(size_t num_threads = 0);
    ~ThreadPool();
    size_t NumThreads() const;
    void Cancel();
    bool IsCancelled() const;
};

}  // namespace OECluster

// ============================================================================
// PDistOptions and pdist function
// Redeclare PDistOptions to handle the std::function member.
// ============================================================================
namespace OECluster {

class DistanceMetric;
class StorageBackend;

struct PDistOptions {
    size_t num_threads;
    size_t chunk_size;
    double cutoff;
    std::function<void(size_t, size_t)> progress;

    PDistOptions();
};

void pdist(DistanceMetric& metric, StorageBackend& storage,
           const PDistOptions& options = PDistOptions());

struct CDistOptions {
    size_t num_threads;
    size_t chunk_size;
    double cutoff;
    std::function<void(size_t, size_t)> progress;

    CDistOptions();
};

void cdist(DistanceMetric& metric, size_t n_a, double* output,
           const CDistOptions& options = CDistOptions());

}  // namespace OECluster

// ============================================================================
// DistanceMatrix
// Redeclare to avoid unique_ptr constructor issues.
// ============================================================================
namespace OECluster {

class DistanceMatrix {
public:
    StorageBackend& Storage();
    const std::string& MetricName() const;
    const std::vector<std::string>& Labels() const;
    size_t NumItems() const;
    size_t NumPairs() const;
};

}  // namespace OECluster

// ============================================================================
// Metric classes
// ============================================================================
%include "oecluster/metrics/FingerprintMetric.h"
%include "oecluster/metrics/ROCSMetric.h"
%include "oecluster/metrics/SuperposeMetric.h"

// ============================================================================
// Version macros
// ============================================================================
#define OECLUSTER_VERSION_MAJOR 3
#define OECLUSTER_VERSION_MINOR 2
#define OECLUSTER_VERSION_PATCH 0

// ============================================================================
// Module-level Python convenience code
// ============================================================================
%pythoncode %{
__version__ = "3.1.7"
%}
