/* -*- mode: c++ -*- */
// swig/oecluster.i
// SWIG interface file for oecluster Python bindings
%module _oecluster

// ============================================================================
// C++ header includes for wrapped API
// ============================================================================
%{
#include "oecluster/oecluster.h"
#include "oecluster/Error.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/GateFacts.h"
#include "oecluster/DescriptorStatistics.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/ThreadPool.h"
#include "oecluster/PDist.h"
#include "oecluster/CDist.h"
#include "oecluster/DistanceMatrix.h"
#include "oecluster/Subset.h"
#include "oecluster/Consensus.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include "oecluster/comparisons/ROCSComparison.h"
#include "oecluster/comparisons/SuperposeComparison.h"
#include "oecluster/comparisons/DescriptorComparison.h"
#include "oecluster/comparisons/RMSDComparison.h"
#include "oecluster/comparisons/MCSComparison.h"
#include "oecluster/clustering/DBSCAN.h"
#include "oecluster/clustering/HDBSCAN.h"
#include "oecluster/clustering/Agglomerative.h"
#include "oecluster/clustering/BitBirch.h"
#include "oecluster/clustering/ISimReport.h"
#include "oecluster/clustering/Representative.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/PartitionAgreement.h"
#include "oecluster/clustering/SARCoherence.h"
#include "oecluster/clustering/DiversitySelection.h"
#include "oecluster/clustering/SetDiversity.h"
#include "oecluster/clustering/SphereExclusion.h"
#include "oecluster/clustering/KNNGraph.h"
#include "oecluster/clustering/JarvisPatrick.h"
#include "oecluster/clustering/Leiden.h"
#include "oecluster/clustering/KMedoids.h"
#include "oecluster/clustering/MurckoScaffold.h"
#include "oefp/batch.h"
#include "oefp/oefp.h"

#include <cstdio>
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

// Returns the pointer an OpenEye proxy wraps. NULL with no exception pending
// means it wraps nothing: it has no ``this``, or a null one. NULL with an
// exception pending means reading ``this`` raised something other than
// AttributeError, which the caller propagates rather than reporting a type.
static void* _oecluster_extract_swig_ptr(PyObject* obj) {
    PyObject* thisAttr = PyObject_GetAttrString(obj, "this");
    if (!thisAttr) {
        if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
            PyErr_Clear();
        }
        return NULL;
    }
    void* ptr = ((_SwigPyObjectCompat*)thisAttr)->ptr;
    Py_DECREF(thisAttr);
    return ptr;
}

// Sets TypeError "<noun> list item <index> is a <type name><tail>", with
// <index> also stored on the instance as ``_oecluster_item_index``. Call it
// with no exception pending. It returns with one pending, that TypeError or
// whatever stopped it being built, for the caller to propagate. The stable
// ABI makes PyTypeObject opaque, so the name is read from the type's __name__
// rather than tp_name.
//
// The attribute marks the error as a conversion's own. cdist renames only a
// marked error, to name the input set (items_named_by_set in
// python/oecluster/_comparisons.py), so an element whose check raised a
// TypeError worded like this one keeps its identity. RMSD's Python item check
// raises a marked copy of this message. A change to the format or the
// attribute must change that module too.
static void _oecluster_list_item_error(const char* noun, Py_ssize_t index, PyObject* item,
                                       const char* tail) {
    // A metaclass can make __name__ raise or return a non-string. An Exception
    // from the read, or a non-string, gives "<unnamed type>": the caller
    // still owes a TypeError naming the position, not whatever reading the
    // name raised. Any other BaseException, such as KeyboardInterrupt or
    // SystemExit, propagates unchanged, because a name lookup must not
    // swallow an interrupt or an exit request. RMSD's Python item check
    // applies the same rule.
    PyObject* name = PyObject_GetAttrString((PyObject*)Py_TYPE(item), "__name__");
    if (!name && !PyErr_ExceptionMatches(PyExc_Exception)) {
        return;
    }
    if (!name || !PyUnicode_Check(name)) {
        PyErr_Clear();
        Py_XDECREF(name);
        name = PyUnicode_FromString("<unnamed type>");
        if (!name) {
            return;
        }
    }
    PyObject* message = PyUnicode_FromFormat("%s list item %zd is a %U%s", noun, index, name,
                                             tail);
    Py_DECREF(name);
    if (!message) {
        return;
    }
    // Built as an instance and marked before it is raised, which needs none
    // of the fetch-and-restore calls deprecated since 3.12.
    PyObject* error = PyObject_CallFunctionObjArgs(PyExc_TypeError, message, NULL);
    Py_DECREF(message);
    if (!error) {
        return;
    }
    PyObject* mark = PyLong_FromSsize_t(index);
    if (!mark || PyObject_SetAttrString(error, "_oecluster_item_index", mark) < 0) {
        // The failure to mark is what surfaces, not the unmarked TypeError:
        // cdist would pass that through with an index into items_a + items_b
        // and no set named, which reads as an index into the caller's set.
        Py_XDECREF(mark);
        Py_DECREF(error);
        return;
    }
    Py_DECREF(mark);
    PyErr_SetObject(PyExc_TypeError, error);
    Py_DECREF(error);
}

// The form every molecule list conversion uses.
static void _oecluster_molecule_item_error(Py_ssize_t index, PyObject* item, const char* tail) {
    _oecluster_list_item_error("molecule", index, item, tail);
}

static PyObject* _oecluster_oefp_batch_type = NULL;
static PyObject* _oecluster_oefp_native_batch_type = NULL;

// Loads and caches oefp's two batch types. Returns 1, or -1 with the import or
// attribute error pending: a broken oefp is a failure to report, not a reason
// to call an object something other than a batch. Both are cached together,
// so a partial failure leaves neither set.
static int _oecluster_load_oefp_batch_types() {
    if (_oecluster_oefp_batch_type && _oecluster_oefp_native_batch_type) {
        return 1;
    }

    PyObject* oefp_mod = PyImport_ImportModule("oefp");
    if (!oefp_mod) {
        return -1;
    }
    PyObject* batch_type = PyObject_GetAttrString(oefp_mod, "OEFPBatch");
    Py_DECREF(oefp_mod);
    if (!batch_type) {
        return -1;
    }

    PyObject* native_mod = PyImport_ImportModule("oefp._native");
    if (!native_mod) {
        Py_DECREF(batch_type);
        return -1;
    }
    PyObject* native_type = PyObject_GetAttrString(native_mod, "_NativeOEFPBatch");
    Py_DECREF(native_mod);
    if (!native_type) {
        Py_DECREF(batch_type);
        return -1;
    }

    _oecluster_oefp_batch_type = batch_type;
    _oecluster_oefp_native_batch_type = native_type;
    return 1;
}

// The same two forms, the same rule for choosing between them, and the same
// limit on what the repeated check can report, as the checkers
// DEFINE_OE_TYPE_CHECKER generates below. An overloaded function that takes
// a batch, such as isim, reads the batch argument in a typecheck first, so it
// is subject to that limit.
static int _oecluster_isinstance_oefp_batch(PyObject* obj) {
    if (_oecluster_load_oefp_batch_types() < 0) {
        return -1;
    }
    int is_batch = PyObject_IsInstance(obj, _oecluster_oefp_batch_type);
    if (is_batch != 0) {
        return is_batch;
    }
    return PyObject_IsInstance(obj, _oecluster_oefp_native_batch_type);
}

static bool _oecluster_may_be_oefp_batch(PyObject* obj) {
    int is_batch = _oecluster_isinstance_oefp_batch(obj);
    if (is_batch < 0) {
        PyErr_Clear();
        return true;
    }
    return is_batch == 1;
}

// Same contract as _oecluster_extract_swig_ptr. A high-level OEFPBatch wraps
// the native batch in ``_native``; a native batch is read directly.
static void* _oecluster_extract_oefp_batch_ptr(PyObject* obj) {
    if (_oecluster_load_oefp_batch_types() < 0) {
        return NULL;
    }
    int is_high_level = PyObject_IsInstance(obj, _oecluster_oefp_batch_type);
    if (is_high_level < 0) {
        return NULL;
    }
    if (!is_high_level) {
        return _oecluster_extract_swig_ptr(obj);
    }
    PyObject* native_obj = PyObject_GetAttrString(obj, "_native");
    if (!native_obj) {
        if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
            PyErr_Clear();
        }
        return NULL;
    }
    void* ptr = _oecluster_extract_swig_ptr(native_obj);
    Py_DECREF(native_obj);
    return ptr;
}

// The installed oefp wheel and this extension each carry their own compiled
// copy of OEFP, and the OEFPBatch typemap hands a raw pointer from one to the
// other. BitBirch then reads that object's words directly, so a layout
// disagreement corrupts memory rather than raising. Compare the version this
// extension compiled against with the version the installed wheel compiled
// against, once, before the first pointer crosses.
static int _oecluster_oefp_abi_state = 0;  // 0 unchecked, 1 ok, -1 mismatch
// Wide enough for the whole diagnostic even when every version component is a
// full-width long; 256 truncated it mid-sentence for ordinary version numbers.
static char _oecluster_oefp_abi_message[512];

static bool _oecluster_check_oefp_abi() {
    if (_oecluster_oefp_abi_state == 1) {
        return true;
    }
    if (_oecluster_oefp_abi_state == -1) {
        PyErr_SetString(PyExc_ImportError, _oecluster_oefp_abi_message);
        return false;
    }

    PyObject* native_mod = PyImport_ImportModule("oefp._native");
    if (!native_mod) {
        return false;  // propagate the import error unchanged
    }

    long installed[3] = {-1, -1, -1};
    const char* names[3] = {"OEFP_VERSION_MAJOR", "OEFP_VERSION_MINOR",
                            "OEFP_VERSION_PATCH"};
    for (int i = 0; i < 3; ++i) {
        PyObject* value = PyObject_GetAttrString(native_mod, names[i]);
        if (!value) {
            Py_DECREF(native_mod);
            return false;
        }
        installed[i] = PyLong_AsLong(value);
        Py_DECREF(value);
        if (installed[i] == -1 && PyErr_Occurred()) {
            Py_DECREF(native_mod);
            return false;
        }
    }
    Py_DECREF(native_mod);

    if (installed[0] == OEFP_VERSION_MAJOR && installed[1] == OEFP_VERSION_MINOR &&
        installed[2] == OEFP_VERSION_PATCH) {
        _oecluster_oefp_abi_state = 1;
        return true;
    }

    std::snprintf(_oecluster_oefp_abi_message, sizeof(_oecluster_oefp_abi_message),
                  "oecluster was compiled against OEFP %d.%d.%d but the installed oefp "
                  "package was compiled against %ld.%ld.%ld. oecluster passes raw "
                  "fingerprint-batch pointers to OEFP, so the two must match exactly. "
                  "Install oefp %d.%d.%d, or rebuild oecluster against the installed "
                  "version.",
                  OEFP_VERSION_MAJOR, OEFP_VERSION_MINOR, OEFP_VERSION_PATCH,
                  installed[0], installed[1], installed[2],
                  OEFP_VERSION_MAJOR, OEFP_VERSION_MINOR, OEFP_VERSION_PATCH);
    _oecluster_oefp_abi_state = -1;
    PyErr_SetString(PyExc_ImportError, _oecluster_oefp_abi_message);
    return false;
}

// ---- Type checker generator macro ----
// Generates a cached isinstance check for an OpenEye Python type, in two forms.
// _oecluster_isinstance_TAG returns 1 or 0, or -1 with the exception set when
// the import or the check itself raised. _oecluster_may_be_TAG clears that
// exception and answers yes, because a check that raised has not ruled the
// object out.
//
// The rule for choosing between them: every %typemap(in) calls
// _oecluster_isinstance_TAG and propagates -1 with SWIG_fail, because a check
// that raised has not answered and its exception is the report. Only a
// %typemap(typecheck) calls _oecluster_may_be_TAG. A typecheck cannot report
// an error, but it can leave its overload in play: if SWIG selects that
// overload, its in typemap repeats the check, and the check's own exception
// reaches the caller only if the repeat raises too. Propagation is therefore
// guaranteed only for a check that answers the same way twice. An object
// whose __class__ or isinstance behaviour changes between consecutive reads
// can raise on the typecheck's read and answer on the in typemap's. Its first
// exception is then lost, and the in typemap goes by the second answer,
// treating the object as that type or refusing it as the wrong one. Only a
// contrived object does that, and the loss is accepted rather than guarded
// against. Answering no would rule the overload out, and when no other
// accepts the argument, replace the exception with the dispatcher's generic
// TypeError.
// TAG:    identifier suffix (e.g., oemolbase)
// MODULE: Python module string (e.g., "openeye.oechem")
// CLASS:  Python class name string (e.g., "OEMolBase")
#define DEFINE_OE_TYPE_CHECKER(TAG, MODULE, CLASS) \
    static PyObject* _oecluster_oe_##TAG##_type = NULL; \
    static int _oecluster_isinstance_##TAG(PyObject* obj) { \
        if (!_oecluster_oe_##TAG##_type) { \
            PyObject* mod = PyImport_ImportModule(MODULE); \
            if (!mod) return -1; \
            _oecluster_oe_##TAG##_type = PyObject_GetAttrString(mod, CLASS); \
            Py_DECREF(mod); \
            if (!_oecluster_oe_##TAG##_type) return -1; \
        } \
        return PyObject_IsInstance(obj, _oecluster_oe_##TAG##_type); \
    } \
    static bool _oecluster_may_be_##TAG(PyObject* obj) { \
        int is_instance = _oecluster_isinstance_##TAG(obj); \
        if (is_instance < 0) { \
            PyErr_Clear(); \
            return true; \
        } \
        return is_instance == 1; \
    }

// ---- Type checker instantiations ----
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
// TAG:      DEFINE_OE_TYPE_CHECKER suffix (e.g., oemolbase); the in typemaps
//           call _oecluster_isinstance_TAG and the typechecks
//           _oecluster_may_be_TAG
// ERR_MSG:  error message on type mismatch
%define OE_CROSS_RUNTIME_REF_TYPEMAPS(CPP_TYPE, TAG, ERR_MSG)

%typemap(in) const CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        int is_instance = _oecluster_isinstance_##TAG($input);
        if (is_instance < 0) {
            SWIG_fail;
        }
        if (is_instance) {
            argp = _oecluster_extract_swig_ptr($input);
            if (argp) {
                res = SWIG_OK;
            } else if (PyErr_Occurred()) {
                SWIG_fail;
            }
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
    $1 = SWIG_IsOK(res) ? 1 : _oecluster_may_be_##TAG($input) ? 1 : 0;
}

%typemap(in) CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        int is_instance = _oecluster_isinstance_##TAG($input);
        if (is_instance < 0) {
            SWIG_fail;
        }
        if (is_instance) {
            argp = _oecluster_extract_swig_ptr($input);
            if (argp) {
                res = SWIG_OK;
            } else if (PyErr_Occurred()) {
                SWIG_fail;
            }
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
    $1 = SWIG_IsOK(res) ? 1 : _oecluster_may_be_##TAG($input) ? 1 : 0;
}

%enddef

// Generate nullable-pointer typemaps (accepts None) for a cross-runtime type.
%define OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(CPP_TYPE, TAG, ERR_MSG)

%typemap(in) const CPP_TYPE* (void *argp = 0, int res = 0) {
    if ($input == Py_None) {
        $1 = NULL;
    } else {
        res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
        if (!SWIG_IsOK(res)) {
            int is_instance = _oecluster_isinstance_##TAG($input);
            if (is_instance < 0) {
                SWIG_fail;
            }
            if (is_instance) {
                argp = _oecluster_extract_swig_ptr($input);
                if (argp) {
                    res = SWIG_OK;
                } else if (PyErr_Occurred()) {
                    SWIG_fail;
                }
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
        $1 = SWIG_IsOK(res) ? 1 : _oecluster_may_be_##TAG($input) ? 1 : 0;
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
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolBase,    oemolbase,    "Expected OEMolBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMCMolBase,  oemcmolbase,  "Expected OEMCMolBase-derived object (OEMCMolBase or OEMol).")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMol,        oemol,        "Expected OEMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEGraphMol,   oegraphmol,   "Expected OEGraphMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEQMol,       oeqmol,       "Expected OEQMol object.")

// ---- Atom / bond / conformer / residue / match (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEAtomBase,   oeatombase,   "Expected OEAtomBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEBondBase,   oebondbase,   "Expected OEBondBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEConfBase,   oeconfbase,   "Expected OEConfBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEResidue,    oeresidue,    "Expected OEResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMatchBase,  oematchbase,  "Expected OEMatchBase-derived object.")

// ---- Molecule I/O (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolistream,  oemolistream, "Expected oemolistream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolostream,  oemolostream, "Expected oemolostream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolDatabase, oemoldatabase,"Expected OEMolDatabase object.")

// ---- Reactions (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEUniMolecularRxn, oeunimolecularrxn, "Expected OEUniMolecularRxn object.")

// ---- Platform streams (OEPlatform) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeifstream, oeifstream, "Expected oeifstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeofstream, oeofstream, "Expected oeofstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeisstream, oeisstream, "Expected oeisstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeosstream, oeosstream, "Expected oeosstream object.")

// ---- Records (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OERecord,    oerecord,    "Expected OERecord object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OEMolRecord, oemolrecord, "Expected OEMolRecord object.")

// ---- Bio / hierarchy (OEBio) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEDesignUnit,   oedesignunit, "Expected OEDesignUnit object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierView,     oehierview,   "Expected OEHierView object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierResidue,  oehierresidue,"Expected OEHierResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierFragment,  oehierfragment,"Expected OEHierFragment object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierChain,    oehierchain,  "Expected OEHierChain object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHint,          oeinteractionhint,          "Expected OEInteractionHint object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHintContainer, oeinteractionhintcontainer, "Expected OEInteractionHintContainer object.")

// ---- Grid (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OEScalarGrid, oescalargrid, "Expected OEScalarGrid-derived object.")
OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(OESystem::OEScalarGrid, oescalargrid, "Expected OEScalarGrid or None.")

// OEScalarGrid return-type typemap (wraps C++ grid as native openeye.oegrid object)
%typemap(out) OESystem::OEScalarGrid* {
    $result = _oecluster_wrap_as_oe_grid($1);
    if (!$result) SWIG_fail;
}

// ---- Docking (OEDocking) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEDocking::OEReceptor, oereceptor, "Expected OEReceptor object.")

// ============================================================================
// Domain-specific typemaps for oecluster's container-based APIs
// ============================================================================
// Molecule and design-unit list conversions. A vector of OEMolBase pointers
// borrows the Python-owned molecules for the call; a vector of shared_ptr
// copies each one, because we do not own the Python-owned objects.
// ROCSComparison and MCSComparison take their lists through a conversion of
// their own, defined after the strict OEMol list conversion.

// The catch clauses for an allocation inside a typemap body. A typemap runs
// outside every %exception handler -- %exception wraps only $action -- so a
// std::bad_alloc escaping one would leave the wrapper as a C++ exception and
// terminate the interpreter. Follow a ``try`` with this. std::bad_alloc and
// std::length_error become SWIG_MemoryError, as in
// OECLUSTER_GIL_MEMORY_EXCEPTION and the consensus out typemaps below, and
// SWIG_exception_fail reaches ``fail:``, so the wrapper's freearg typemaps
// still run. WHAT completes "out of memory ...".
%define OECLUSTER_TYPEMAP_ALLOCATION_CATCH(WHAT)
catch (const std::bad_alloc&) {
    SWIG_exception_fail(SWIG_MemoryError, "out of memory " WHAT);
} catch (const std::length_error& e) {
    SWIG_exception_fail(SWIG_MemoryError, e.what());
} catch (const std::exception& e) {
    SWIG_exception_fail(SWIG_RuntimeError, e.what());
} catch (...) {
    SWIG_exception_fail(SWIG_RuntimeError, "Unknown C++ exception " WHAT);
}
%enddef

// ---- const std::vector<OEChem::OEMolBase*>& ----
%typemap(in) const std::vector<OEChem::OEMolBase*>& (std::vector<OEChem::OEMolBase*> temp, PyObject* held = 0) {
    if (!PyList_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError, "Expected a list of OEMolBase objects");
    }
    // The tuple holds a reference to every element until the freearg typemap
    // below releases it, so for the duration of the call, emptying the
    // caller's list cannot free a molecule. An entry point that releases the
    // GIL lets another thread do exactly that mid-call. Nothing here protects
    // a pointer kept after the call returns. The loop reads the tuple rather
    // than the list because its checks can run Python code, which could
    // resize the list under it.
    held = PyList_AsTuple($input);
    if (!held) {
        SWIG_fail;
    }
    Py_ssize_t size = PyTuple_GET_SIZE(held);
    try {
        temp.reserve(size);
    } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a molecule list")
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyTuple_GET_ITEM(held, i);
        // A check that raised has not answered, so its own exception is the
        // report, not a claim about the item's type.
        int is_mol = _oecluster_isinstance_oemolbase(item);
        if (is_mol < 0) {
            SWIG_fail;
        }
        if (!is_mol) {
            _oecluster_molecule_item_error(i, item, ", not an OEMolBase");
            SWIG_fail;
        }
        // An OEMolBase instance whose ``this`` is missing or null, such as one
        // made by ``__new__`` without ``__init__``, passes the check above, so
        // it gets its own message rather than being called a non-OEMolBase.
        void* ptr = _oecluster_extract_swig_ptr(item);
        if (!ptr) {
            if (!PyErr_Occurred()) {
                _oecluster_molecule_item_error(i, item, " that wraps no molecule");
            }
            SWIG_fail;
        }
        // Same offset-0 base-cast assumption as the shared_ptr<OEMol> typemap
        // below; see the comment there. The reserve above means this cannot
        // reallocate today; the guard keeps that from being load-bearing.
        try {
            temp.push_back(reinterpret_cast<OEChem::OEMolBase*>(ptr));
        } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a molecule list")
    }
    $1 = &temp;
}

// Py_XDECREF needs the GIL. This relies on every %exception handler having
// reacquired it before the wrapper reaches freearg, on the fail: path as well
// as on success.
%typemap(freearg) const std::vector<OEChem::OEMolBase*>& {
    Py_XDECREF(held$argnum);
}

// ``precedence=1``, a literal rather than ``SWIG_TYPECHECK_POINTER``, and
// deliberate. A bare numeric is ordinary here -- four typechecks in this file
// use ``precedence=10`` -- so the value is not an oddity to be tidied into a
// symbolic constant.
//
// What it orders is SuperposeComparison's two constructors: one takes this
// vector, the other the design-unit vector below, whose typecheck sits at
// ``SWIG_TYPECHECK_POINTER``, which is 0. SWIG ranks overloads by argument
// count first and typemap precedence second, and lower precedence is examined
// first, so the design-unit typecheck is examined ahead of this one in every
// argument-count group. Equal precedence would hand the order to SWIG's
// unspecified tie-break, which was measured, on the ROCS constructors this
// value used to order, to rank a pair one way in the one-argument group and
// the opposite way in the two-argument group.
//
// This typecheck tests list-ness only, which makes SuperposeComparison's
// correctness depend on that order. FingerprintComparison,
// DescriptorComparison and SuperposeComparison each keep one- and
// two-argument constructor forms, so SWIG emits a dispatcher in front of
// them, and a typecheck that rejected a bad element there would replace the
// conversion's positional message with the dispatcher's generic
// "new_<Class>" one, which names no element. (A function that takes this
// vector with no overload has no dispatcher, so the typecheck is never
// consulted for it.) Accepting every list means that, were this typecheck
// examined first, it would take design-unit lists too, and the conversion
// above would refuse each design unit as "not an OEMolBase". The
// design-unit typecheck still reads element 0 (its comment says what that
// read can do to the list), and being examined first is what keeps a list
// that starts with a design unit on the design-unit constructor, leaving
// every other list to this one. Two lists reach the
// design-unit constructor without element 0 being known to be a design unit:
// an empty list, refused with the same "SuperposeComparison: empty structure
// list" either constructor gives; and one whose first element's isinstance
// check raises, which the design-unit typecheck leaves in play, so the
// design-unit conversion is the one that repeats the check. Fingerprint and
// Descriptor have no competing overload for the order to matter to.
//
// ROCSComparison and MCSComparison no longer take part. Each used to carry a
// constructor taking this vector beside one taking the strict
// ``shared_ptr<OEMol>`` vector, and this value was what kept a list of OEMol on
// the conformer-preserving one. Their constructors taking this vector are now
// %ignore'd, and both classes take every molecule list through the one
// conversion defined after the strict typemap below.
%typemap(typecheck, precedence=1) const std::vector<OEChem::OEMolBase*>& {
    $1 = PyList_Check($input) ? 1 : 0;
}

// ---- const std::vector<std::shared_ptr<OEChem::OEMol>>& ----
%typemap(in) const std::vector<std::shared_ptr<OEChem::OEMol>>&
        (std::vector<std::shared_ptr<OEChem::OEMol>> temp, PyObject* held = 0) {
    if (!PyList_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError, "Expected a list of OEMol objects");
    }
    // Held, and read, for the same reasons as in the OEMolBase list
    // conversion above, and released by the freearg below. The message for a
    // refused element reads that element's type after its check has run, and
    // the check can run Python code that empties the caller's list.
    held = PyList_AsTuple($input);
    if (!held) {
        SWIG_fail;
    }
    Py_ssize_t size = PyTuple_GET_SIZE(held);
    try {
        temp.reserve(size);
    } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a molecule list")
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyTuple_GET_ITEM(held, i);
        int is_instance = _oecluster_isinstance_oemol(item);
        if (is_instance < 0) {
            SWIG_fail;
        }
        if (!is_instance) {
            _oecluster_molecule_item_error(i, item, ", not an OEMol");
            SWIG_fail;
        }
        void* ptr = _oecluster_extract_swig_ptr(item);
        if (!ptr) {
            if (!PyErr_Occurred()) {
                _oecluster_molecule_item_error(i, item, " that wraps no molecule");
            }
            SWIG_fail;
        }
        // Make a COPY wrapped in shared_ptr (we do not own the Python object).
        //
        // The cast target is OEMCMolBase, not OEMol. OpenEye's Python bindings
        // hand out an OEMolWrapper*, which derives from OEMCMolBase but is not
        // an OEChem::OEMol -- that type is a bare handle holding an
        // OEOwnedPtr<OEMCMolBase>. Casting to OEMol and copy-constructing reads
        // the wrapper's first member as if it were that owned pointer and
        // segfaults. OEMol(const OEMCMolBase&) is the documented converting
        // constructor for exactly this, and it preserves conformers.
        //
        // The cast is a base-class cast through a type-erased void*, so it is
        // correct only while the OEMCMolBase subobject sits at offset 0 in
        // OEMolWrapper. That holds today -- OEMCMolBase -> OEMolBase ->
        // OESystem::OEBase is single, non-virtual inheritance throughout -- but
        // multiple or virtual bases would shift the subobject, and this cast
        // would then read the wrong bytes without crashing. SWIG's own base
        // conversion cannot be used instead: the pointer comes from another
        // extension module's type registry, which this module cannot query.
        auto* mcmol = reinterpret_cast<OEChem::OEMCMolBase*>(ptr);
        try {
            temp.push_back(std::make_shared<OEChem::OEMol>(*mcmol));
        } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a molecule list")
    }
    $1 = &temp;
}

// Same GIL precondition as the OEMolBase list's freearg.
%typemap(freearg) const std::vector<std::shared_ptr<OEChem::OEMol>>& {
    Py_XDECREF(held$argnum);
}

// List-ness only, for the reason given at the OEMolBase list's typecheck:
// RMSDComparison has one- and two-argument constructor forms, and a
// typecheck that rejected a bad element would replace the positional message
// above with the dispatcher's.
%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) const std::vector<std::shared_ptr<OEChem::OEMol>>& {
    $1 = PyList_Check($input) ? 1 : 0;
}

// ---- ROCSComparison and MCSComparison: any OEMolBase, each by its own type ----
// Each element is copied by its own type: an OEMCMolBase, which is what an
// OEMol is, through ``OEMol(const OEMCMolBase&)``, keeping every conformer;
// anything else through ``OEMol(const OEMolBase&)``, as the one conformer it
// has. No element's type decides how another is read, so an element is read
// the same wherever it sits. That used not to hold: each class also exposed a
// constructor taking ``vector<OEMolBase*>``, which copies through the second
// constructor alone, and the typechecks chose between the two by element 0,
// so a list starting with an OEGraphMol scored every later OEMol on its
// active conformer.
//
// RMSDComparison takes the same parameter type and keeps the strict
// conversion above, which refuses an OEGraphMol. That is a scope decision,
// not a scoring requirement: RMSD reads only each molecule's active
// conformer, which an OEGraphMol supplies, but widening its conversion would
// change what it accepts, and that is a separate decision. So this is defined
// under a parameter name nothing declares and %apply'd to
// ``mols`` only around the ROCSComparison.h and MCSComparison.h %includes,
// %clear'ed after each. RMSDComparison.h is %included between those two,
// which is why there are two windows rather than one. The windows key on the
// parameter name: renaming ``mols`` in either header puts that class back on
// the strict conversion with a clean build, which only the Python tests would
// notice.
//
// The typecheck tests list-ness only. Both classes keep one- and
// two-argument constructor forms, so SWIG still emits a dispatcher, and a
// typecheck that rejected a bad element there would replace this conversion's
// positional message with the dispatcher's generic "new_<Class>" one.
%typemap(in) const std::vector<std::shared_ptr<OEChem::OEMol>>& ANY_MOLECULE_LIST
        (std::vector<std::shared_ptr<OEChem::OEMol>> temp, PyObject* held = 0) {
    if (!PyList_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError, "Expected a list of OEMolBase objects");
    }
    // Held, and read, for the same reasons as in the OEMolBase list
    // conversion above, and released by the freearg below. Every element is
    // copied before the call today; holding the tuple for the whole call
    // anyway means a change that borrows instead of copying cannot reopen the
    // lifetime hole that conversion closes.
    held = PyList_AsTuple($input);
    if (!held) {
        SWIG_fail;
    }
    Py_ssize_t size = PyTuple_GET_SIZE(held);
    try {
        temp.reserve(size);
    } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a molecule list")
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyTuple_GET_ITEM(held, i);
        int is_mol = _oecluster_isinstance_oemolbase(item);
        if (is_mol < 0) {
            SWIG_fail;
        }
        if (!is_mol) {
            _oecluster_molecule_item_error(i, item, ", not an OEMolBase");
            SWIG_fail;
        }
        int is_multi = _oecluster_isinstance_oemcmolbase(item);
        if (is_multi < 0) {
            SWIG_fail;
        }
        void* ptr = _oecluster_extract_swig_ptr(item);
        if (!ptr) {
            if (!PyErr_Occurred()) {
                _oecluster_molecule_item_error(i, item, " that wraps no molecule");
            }
            SWIG_fail;
        }
        // The same offset-0 base casts as the two conversions above. The copy
        // runs before $action, outside every %exception handler, so it is
        // guarded here by OECLUSTER_TYPEMAP_ALLOCATION_CATCH.
        try {
            if (is_multi) {
                temp.push_back(std::make_shared<OEChem::OEMol>(
                    *reinterpret_cast<OEChem::OEMCMolBase*>(ptr)));
            } else {
                temp.push_back(std::make_shared<OEChem::OEMol>(
                    *reinterpret_cast<OEChem::OEMolBase*>(ptr)));
            }
        } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a molecule list")
    }
    $1 = &temp;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) const std::vector<std::shared_ptr<OEChem::OEMol>>& ANY_MOLECULE_LIST {
    $1 = PyList_Check($input) ? 1 : 0;
}

// Same GIL precondition as the OEMolBase list's freearg.
%typemap(freearg) const std::vector<std::shared_ptr<OEChem::OEMol>>& ANY_MOLECULE_LIST {
    Py_XDECREF(held$argnum);
}

// ---- const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& ----
%typemap(in) const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>&
        (std::vector<std::shared_ptr<OEBio::OEDesignUnit>> temp, PyObject* held = 0) {
    if (!PyList_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError, "Expected a list of OEDesignUnit objects");
    }
    // Held, read and released as in the strict OEMol list conversion above.
    held = PyList_AsTuple($input);
    if (!held) {
        SWIG_fail;
    }
    Py_ssize_t size = PyTuple_GET_SIZE(held);
    try {
        temp.reserve(size);
    } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a design-unit list")
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyTuple_GET_ITEM(held, i);
        int is_instance = _oecluster_isinstance_oedesignunit(item);
        if (is_instance < 0) {
            SWIG_fail;
        }
        if (!is_instance) {
            _oecluster_list_item_error("design unit", i, item, ", not an OEDesignUnit");
            SWIG_fail;
        }
        void* ptr = _oecluster_extract_swig_ptr(item);
        if (!ptr) {
            if (!PyErr_Occurred()) {
                _oecluster_list_item_error("design unit", i, item, " that wraps no design unit");
            }
            SWIG_fail;
        }
        // Make a COPY wrapped in shared_ptr (we do not own the Python object)
        auto* du = reinterpret_cast<OEBio::OEDesignUnit*>(ptr);
        try {
            temp.push_back(std::make_shared<OEBio::OEDesignUnit>(*du));
        } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a design-unit list")
    }
    $1 = &temp;
}

// Same GIL precondition as the OEMolBase list's freearg.
%typemap(freearg) const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& {
    Py_XDECREF(held$argnum);
}

// Unlike the molecule list typechecks, this one still reads element 0: it is
// what chooses between SuperposeComparison's two constructors. See the
// comment at the OEMolBase list's typecheck for the order that relies on.
//
// What the caller's list is, for every molecule and design-unit list
// conversion in this file: overload resolution can run Python code before
// any conversion does. This check runs element 0's isinstance, which can
// read its __class__, and an options argument's own type check can read its
// ``this``. Either can change the list. The conversion SWIG then selects
// snapshots and checks the list as it stands when that conversion runs, not
// as it stood at call entry, and every element the conversion reads is
// still checked and held for the call.
//
// Element 0 is held across the check because a __class__ that empties the
// list can drop its last reference while isinstance is still reading it.
%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& {
    if (PyList_Check($input) && PyList_Size($input) > 0) {
        PyObject* first = PyList_GetItem($input, 0);
        Py_XINCREF(first);
        $1 = _oecluster_may_be_oedesignunit(first) ? 1 : 0;
        Py_XDECREF(first);
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
        // Both statements below can allocate: the shared_ptr its control
        // block, and the std::function its copy of the lambda when that is
        // stored on the heap. Neither leaks the reference taken above if it
        // throws. shared_ptr calls the deleter itself when the control block
        // cannot be allocated, and a throwing assignment unwinds through
        // guard's destructor, before the catch runs.
        try {
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
        } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a progress callback")
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

// ---- Return by value ----
// SWIG's stock typemap copies the returned value onto the heap after
// %exception has closed, so a bad_alloc there would escape into the
// interpreter and abort it. Moving instead of copying also avoids holding two
// full copies of a large result at once. This must precede the STL includes,
// whose container templates are the first by-value returns.
%typemap(out, noblock=1) SWIGTYPE {
  try {
    $result = SWIG_NewPointerObj((new $1_ltype(SWIG_STD_MOVE($1))), $&descriptor,
                                 SWIG_POINTER_OWN | %newpointer_flags);
  } catch (const std::bad_alloc&) {
    SWIG_exception_fail(SWIG_MemoryError, "out of memory returning $1_basetype");
  } catch (const std::exception& e) {
    SWIG_exception_fail(SWIG_RuntimeError, e.what());
  } catch (...) {
    SWIG_exception_fail(SWIG_RuntimeError, "Unknown C++ exception returning $1_basetype");
  }
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
%template(SizeTVector) std::vector<size_t>;
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(ClusterVector) std::vector<std::vector<size_t>>;

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
// GIL release for the long-running entry points (overrides the general
// handler above)
//
// Every override releases the GIL around $action and restores it on each exit
// path. Most of them then share one catch ladder, differing only in the name
// in the fallback message; spelling that out by hand for each one invited the
// ladders to drift apart, and one of them already had. `NAME` is SWIG
// preprocessor stringification, which is what keeps each message naming its
// own function. The entry points that can plausibly exhaust memory use
// OECLUSTER_GIL_MEMORY_EXCEPTION instead, which adds a pair of catch clauses.
// ============================================================================
%define OECLUSTER_GIL_EXCEPTION(QUALIFIED, NAME)
%exception QUALIFIED {
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
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception in " `NAME`);
    }
    Py_END_ALLOW_THREADS
}
%enddef

// The same ladder, except that std::bad_alloc and std::length_error (a
// container asked to exceed its max_size()) surface as SWIG_MemoryError rather
// than the blanket SWIG_RuntimeError, so a caller can tell an exhausted
// machine from a refused input.
%define OECLUSTER_GIL_MEMORY_EXCEPTION(QUALIFIED, NAME)
%exception QUALIFIED {
    Py_BEGIN_ALLOW_THREADS
    try {
        $action
    } catch (const OECluster::OEClusterError& e) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (const std::bad_alloc&) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_MemoryError, `NAME` ": out of memory");
    } catch (const std::length_error& e) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_MemoryError, e.what());
    } catch (const std::exception& e) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        Py_BLOCK_THREADS
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception in " `NAME`);
    }
    Py_END_ALLOW_THREADS
}
%enddef

OECLUSTER_GIL_EXCEPTION(OECluster::pdist, pdist)
OECLUSTER_GIL_EXCEPTION(OECluster::cdist_into_address, cdist_into_address)
OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::butina_cluster, butina_cluster)
OECLUSTER_GIL_EXCEPTION(OECluster::cluster_representative, cluster_representative)
OECLUSTER_GIL_EXCEPTION(OECluster::rank_representatives, rank_representatives)
OECLUSTER_GIL_EXCEPTION(OECluster::select_representatives, select_representatives)

// cluster_report allocates proportionally to the square of the input (its
// storage overload under compute_pair_rank_indices), so it can plausibly
// exhaust memory or exceed a container's max_size(). The handler covers both
// overloads; the comparison overload runs clones on worker threads with no
// Python object in reach, so the GIL stays released.
OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::cluster_report, cluster_report)

OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::dbscan_cluster, dbscan_cluster)
// Both comparison overloads hold O(N) spanning-tree arrays, and HDBSCAN's
// holds N x (min_samples - 1) core-distance heaps, so an exhausted machine
// reads as MemoryError on every path, the matrix paths included.
OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::hdbscan_cluster, hdbscan_cluster)
OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::agglomerative_cluster, agglomerative_cluster)
OECLUSTER_GIL_EXCEPTION(OECluster::bitbirch_cluster, bitbirch_cluster)
// isim_report's centroid stage holds a K * bits * 4 byte column store, about
// 0.8 GB for 10^5 clusters of 2048 bits, which is within reach of a BitBirch
// clustering. isim shares the handler so the pair refuses alike.
OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::isim, isim)
OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::isim_report, isim_report)
OECLUSTER_GIL_EXCEPTION(OECluster::bitbirch_recluster, bitbirch_recluster)
OECLUSTER_GIL_EXCEPTION(OECluster::bitbirch_refine, bitbirch_refine)

// BUILD initialization at O(k N^2) and the swap loop's per-iteration pass both
// sweep a precomputed matrix, and the ThreadPool workers the call spawns run
// entirely on native data with no Python object in reach.
OECLUSTER_GIL_EXCEPTION(OECluster::k_medoids_cluster, k_medoids_cluster)

// Extraction spawns ThreadPool workers that touch only OEChem data, and
// nothing it builds touches a Python object before the GIL is reacquired, so
// the GIL is released for both entry points.
OECLUSTER_GIL_EXCEPTION(OECluster::murcko_scaffolds, murcko_scaffolds)
OECLUSTER_GIL_EXCEPTION(OECluster::murcko_cluster, murcko_cluster)

// The AMI path is the library's one second-scale computation, so it must not
// hold the GIL. The override is chosen at wrap time and cannot see the option
// value, so it wraps every overload of the name rather than branching.
OECLUSTER_GIL_EXCEPTION(OECluster::partition_agreement, partition_agreement)
OECLUSTER_GIL_EXCEPTION(OECluster::scaffold_agreement, scaffold_agreement)

// activity_landscape and modelability sweep every pair, of a precomputed
// matrix or of a comparison whose clones run on worker threads; sar_coherence
// is linear. None of them has a Python object in reach, so the GIL stays
// released across every overload.
OECLUSTER_GIL_EXCEPTION(OECluster::sar_coherence, sar_coherence)
OECLUSTER_GIL_EXCEPTION(OECluster::activity_landscape, activity_landscape)
OECLUSTER_GIL_EXCEPTION(OECluster::modelability, modelability)

// maxmin_select and circles fold O(N) native rows per pick, and their
// comparison overloads run for as long as the comparisons do.
OECLUSTER_GIL_EXCEPTION(OECluster::maxmin_select, maxmin_select)
OECLUSTER_GIL_EXCEPTION(OECluster::circles, circles)

// vendi_score and logdet_diversity build and decompose an O(N^2) kernel, or
// at order 2 run O(N^2) comparisons, with no Python object in reach.
OECLUSTER_GIL_EXCEPTION(OECluster::vendi_score, vendi_score)
OECLUSTER_GIL_EXCEPTION(OECluster::logdet_diversity, logdet_diversity)

// sphere_exclusion reads O(N*k) native rows, or builds the O(N^2) threshold
// graph under the neighbor order, and its comparison overload runs clones on
// worker threads for as long as the comparisons do.
OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::sphere_exclusion, sphere_exclusion)

// knn_graph reads O(N^2) native distances or runs O(N^2) comparisons on worker
// threads, and jarvis_patrick builds the same graph before linking it.
OECLUSTER_GIL_EXCEPTION(OECluster::knn_graph, knn_graph)
OECLUSTER_GIL_EXCEPTION(OECluster::jarvis_patrick, jarvis_patrick)

// take_pairs reads m(m-1)/2 native distances on worker threads and
// take_fingerprints copies m rows; neither has a Python object in reach, and
// the input typemaps finished with the GIL held before $action runs.
OECLUSTER_GIL_EXCEPTION(OECluster::take_pairs, take_pairs)
OECLUSTER_GIL_EXCEPTION(OECluster::take_fingerprints, take_fingerprints)

// coassociation_distances accumulates R passes over O(N^2) pairs on worker
// threads, consensus_components sweeps the whole condensed array, and
// consensus_strength is another O(N^2) pass; none has a Python object in
// reach once the typemaps have run.
OECLUSTER_GIL_EXCEPTION(OECluster::coassociation_distances, coassociation_distances)
OECLUSTER_GIL_EXCEPTION(OECluster::consensus_components, consensus_components)
OECLUSTER_GIL_EXCEPTION(OECluster::consensus_strength, consensus_strength)

// leiden builds the same graph, weighs it on worker threads and runs a
// serial optimizer that can take seconds at scale.
OECLUSTER_GIL_EXCEPTION(OECluster::leiden, leiden)

// ============================================================================
// Ignore problematic members before %include
//
// StorageBackend.h pulls in <shared_mutex>, <unordered_map>, <thread>
// which SWIG cannot parse. We redeclare the classes manually instead.
//
// Also ignore private clone constructors for comparisons.
// ============================================================================

// Suppress default constructors for types that don't have them
%nodefaultctor OECluster::OEClusterError;
%nodefaultctor OECluster::ComparisonError;
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

// Ignore PairwiseComparison::Clone -- returns unique_ptr, expose via %extend
%ignore OECluster::PairwiseComparison::Clone;
%ignore OECluster::PairwiseComparison::TryPDist;
%ignore OECluster::PairwiseComparison::TryCDist;

// Ignore private constructors for comparison clones
%ignore OECluster::FingerprintComparison::FingerprintComparison(std::shared_ptr<const Impl>);
%ignore OECluster::FingerprintComparison::TryPDist;
%ignore OECluster::FingerprintComparison::TryCDist;

%ignore OECluster::ROCSComparison::ROCSComparison(std::shared_ptr<const SharedData>, const Options&);
%ignore OECluster::SuperposeComparison::SuperposeComparison(std::shared_ptr<const SharedData>, const Options&);

%ignore OECluster::DescriptorComparison::DescriptorComparison(std::shared_ptr<const Impl>);
%ignore OECluster::DescriptorComparison::TryPDist;
%ignore OECluster::DescriptorComparison::TryCDist;

%ignore OECluster::RMSDComparison::RMSDComparison(std::shared_ptr<const SharedData>, const Options&);
%ignore OECluster::MCSComparison::MCSComparison(std::shared_ptr<const SharedData>, const Options&);

// Hidden from Python only; C++ callers keep both. The remaining constructor's
// list conversion (ANY_MOLECULE_LIST) accepts every OEMolBase, so these would
// only offer a second, active-conformer reading of the same list.
%ignore OECluster::ROCSComparison::ROCSComparison(const std::vector<OEChem::OEMolBase*>&, const Options& = Options());
%ignore OECluster::MCSComparison::MCSComparison(const std::vector<OEChem::OEMolBase*>&, const Options& = Options());

// Ignore SparseStorage internals that use unordered_map/shared_mutex/thread
%ignore OECluster::SparseStorage::Entries;

// The Python KNNGraph is built only by knn_graph(), which maps rows to caller
// positions; a graph adopted from raw arrays would carry no such mapping.
%ignore OECluster::KNNGraph::KNNGraph(size_t, size_t, std::vector<size_t>, std::vector<double>);

// Accept high-level oefp.OEFPBatch and native oefp._native._NativeOEFPBatch
// objects without making oecluster own or duplicate fingerprint wrappers.
namespace OEFP {
class OEFPBatch;
}

%typemap(in) const OEFP::OEFPBatch& (void *argp = 0) {
    if (!_oecluster_check_oefp_abi()) {
        SWIG_fail;
    }
    int is_batch = _oecluster_isinstance_oefp_batch($input);
    if (is_batch < 0) {
        SWIG_fail;
    }
    if (!is_batch) {
        SWIG_exception_fail(SWIG_TypeError, "Expected an oefp.OEFPBatch.");
    }
    argp = _oecluster_extract_oefp_batch_ptr($input);
    if (!argp) {
        if (PyErr_Occurred()) {
            SWIG_fail;
        }
        SWIG_exception_fail(SWIG_NullReferenceError, "Null OEFPBatch reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) const OEFP::OEFPBatch& {
    $1 = _oecluster_may_be_oefp_batch($input) ? 1 : 0;
}

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
    virtual size_t NumSamples() const = 0;
    virtual size_t NumPairs() const = 0;
    virtual void Finalize();
};

class DenseStorage : public StorageBackend {
public:
    explicit DenseStorage(size_t n);
    void Set(size_t i, size_t j, double value) override;
    double Get(size_t i, size_t j) const override;
    size_t NumSamples() const override;
    size_t NumPairs() const override;
};

class MMapStorage : public StorageBackend {
public:
    MMapStorage(const std::string& path, size_t n);
    ~MMapStorage() override;
    void Set(size_t i, size_t j, double value) override;
    double Get(size_t i, size_t j) const override;
    size_t NumSamples() const override;
    size_t NumPairs() const override;
};

class SparseStorage : public StorageBackend {
public:
    SparseStorage(size_t n, double cutoff);
    ~SparseStorage() override;
    void Set(size_t i, size_t j, double value) override;
    double Get(size_t i, size_t j) const override;
    size_t NumSamples() const override;
    size_t NumPairs() const override;
    double Cutoff() const;
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
// Capability facts (must precede PairwiseComparison, which returns GateFacts)
// ============================================================================
%include "oecluster/GateFacts.h"

// ============================================================================
// PairwiseComparison base class
// ============================================================================
%include "oecluster/PairwiseComparison.h"

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

class PairwiseComparison;
class StorageBackend;

struct PDistOptions {
    size_t num_threads;
    size_t chunk_size;
    double cutoff;
    std::function<void(size_t, size_t)> progress;

    PDistOptions();
};

void pdist(PairwiseComparison& comparison, StorageBackend& storage,
           const PDistOptions& options = PDistOptions());

struct CDistOptions {
    size_t num_threads;
    size_t chunk_size;
    double cutoff;
    std::function<void(size_t, size_t)> progress;

    CDistOptions();
};

}  // namespace OECluster

%inline %{
namespace OECluster {
// Buffer-address bridge so Python can pass a preallocated numpy array's data
// pointer (as an integer) instead of an opaque double*. The Python wrapper owns
// all dtype/shape/lifetime validation; this only forwards the address.
void cdist_into_address(PairwiseComparison& comparison, size_t n_a,
                        size_t out_addr, const CDistOptions& options) {
    cdist(comparison, n_a, reinterpret_cast<double*>(out_addr), options);
}
}  // namespace OECluster
%}

// ============================================================================
// DistanceMatrix
// Redeclare to avoid unique_ptr constructor issues.
// ============================================================================
namespace OECluster {

class DistanceMatrix {
public:
    StorageBackend& Storage();
    const std::string& ComparisonName() const;
    const std::vector<std::string>& Labels() const;
    size_t NumSamples() const;
    size_t NumPairs() const;
};

}  // namespace OECluster

// ============================================================================
// Comparison classes
// ============================================================================
%include "oecluster/comparisons/FingerprintComparison.h"
// One of the two windows for ANY_MOLECULE_LIST; see its definition for why
// RMSD must stay outside them.
%apply const std::vector<std::shared_ptr<OEChem::OEMol>>& ANY_MOLECULE_LIST {
    const std::vector<std::shared_ptr<OEChem::OEMol>>& mols
};
%include "oecluster/comparisons/ROCSComparison.h"
%clear const std::vector<std::shared_ptr<OEChem::OEMol>>& mols;
%include "oecluster/comparisons/SuperposeComparison.h"
%include "oecluster/comparisons/DescriptorComparison.h"
%include "oecluster/comparisons/RMSDComparison.h"

// ---- MCSOptions::max_matches: refuse a Python bool at the setter ----
// Python's bool is an int subclass, so SWIG's default unsigned int conversion
// accepts True and stores 1. One is the most destructive value this option
// has: a single-match budget changes scores silently instead of raising. The
// refusal therefore belongs at the setter boundary, which every route crosses.
// The wrapper's own guard in _comparisons.mcs_options covers pdist, cdist and
// oecluster.MCSComparison, but a caller who builds the struct through
// oecluster.oecluster -- the access path documented for the option structs
// that are not exported -- and hands the comparison straight to pdist never
// passes through it. That guard still fires first on the wrapper paths, where
// it refuses before any molecule is read and words the error better; this is
// the floor beneath it, not a replacement.
//
// The typemap must be keyed on the bare member name and bracketed by %clear,
// rather than written as `unsigned int OECluster::MCSOptions::max_matches`.
// For a member setter SWIG searches `in` typemaps only under the unqualified
// name -- `swig -debug-tmsearch` reports `Looking for: unsigned int
// max_matches` and nothing else -- so the qualified form matches nothing and
// is accepted in silence, with no warning and a clean build. The qualified
// spelling is consulted only for the getter's `out` typemap.
//
// Hence the placement: the pair brackets exactly one %include, so the typemap
// is live only while MCSComparison.h is parsed. That keeps it to this one
// member deliberately. Sibling option structs stay lenient about bool, which
// is a known inconsistency held for a separate decision, not an oversight.
%typemap(in, fragment="SWIG_AsVal_unsigned_SS_int") unsigned int max_matches {
    if (PyBool_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError,
            "in method '$symname', argument $argnum of type '$type': "
            "max_matches must be an integer, not a bool; True would be read "
            "as a match budget of 1, which silently corrupts scores");
    }
    // Otherwise defer to the conversion SWIG would have generated, so floats,
    // strings, numpy scalars and negatives keep the exact exception types and
    // messages they raised before this typemap existed.
    unsigned int val = 0;
    int ecode = SWIG_AsVal_unsigned_SS_int($input, &val);
    if (!SWIG_IsOK(ecode)) {
        SWIG_exception_fail(SWIG_ArgError(ecode),
            "in method '$symname', argument $argnum of type '$type'");
    }
    $1 = static_cast<unsigned int>(val);
}
// The other ANY_MOLECULE_LIST window.
%apply const std::vector<std::shared_ptr<OEChem::OEMol>>& ANY_MOLECULE_LIST {
    const std::vector<std::shared_ptr<OEChem::OEMol>>& mols
};
%include "oecluster/comparisons/MCSComparison.h"
%clear const std::vector<std::shared_ptr<OEChem::OEMol>>& mols;
%clear unsigned int max_matches;

// ============================================================================
// Descriptor statistics
// ============================================================================
%include "oecluster/DescriptorStatistics.h"

// ============================================================================
// Clustering algorithms
// ============================================================================
%include "oecluster/clustering/ClusterTypes.h"
%include "oecluster/clustering/Butina.h"
%include "oecluster/clustering/Representative.h"
%template(ClusterRepresentativeVector) std::vector<OECluster::ClusterRepresentative>;
%template(ClusterRecordVector) std::vector<OECluster::ClusterRecord>;
%include "oecluster/clustering/ClusterReport.h"
%include "oecluster/clustering/PartitionAgreement.h"
// Declared ahead of their header, as ClusterRecordVector is above: both name
// a vector that appears as a member of a result struct in the header that
// follows. SWIG accepts either order here -- it emits the same wrapper, only
// positioned differently in the generated module -- so this is a convention
// rather than a requirement.
%template(ClusterActivityVector) std::vector<OECluster::ClusterActivity>;
%template(ClassConcordanceVector) std::vector<OECluster::ClassConcordance>;
%include "oecluster/clustering/SARCoherence.h"
%include "oecluster/clustering/DiversitySelection.h"
%include "oecluster/clustering/SetDiversity.h"
%include "oecluster/clustering/SphereExclusion.h"
%include "oecluster/clustering/KNNGraph.h"
%include "oecluster/clustering/JarvisPatrick.h"
%include "oecluster/clustering/Leiden.h"
%include "oecluster/clustering/DBSCAN.h"
%include "oecluster/clustering/HDBSCAN.h"
%include "oecluster/clustering/Agglomerative.h"
%include "oecluster/clustering/BitBirch.h"
// Declared ahead of its header, as ClusterRecordVector is: the vector is a
// member of a result struct the header defines.
%template(ISimClusterRecordVector) std::vector<OECluster::ISimClusterRecord>;
%include "oecluster/clustering/ISimReport.h"
%include "oecluster/clustering/KMedoids.h"
// Both Murcko entry points default their options argument, and SWIG's usual
// expansion of a default argument into two overloads puts an overload
// dispatcher in front of them. When that dispatcher rejects a call it reports
// its own message -- "Wrong number or type of arguments for overloaded function
// 'murcko_cluster'" -- in place of the typemap's, naming a C++ symbol that
// murcko() does not expose. The molecule-list typecheck it consults once read
// the first list element, so a bad first element got that message, which said
// nothing about which element was wrong; it now tests list-ness only, which
// still leaves a non-list to the dispatcher. Compacting the defaults into one
// wrapper each removes the dispatcher, so the molecule-list typemap's own
// message reaches the caller for every input it refuses, as it already does
// for the unoverloaded descriptor_statistics.
%feature("compactdefaultargs") OECluster::murcko_scaffolds;
%feature("compactdefaultargs") OECluster::murcko_cluster;
%include "oecluster/clustering/MurckoScaffold.h"

// ============================================================================
// Row subsets
// ============================================================================
// take_pairs defaults num_threads and chunk_size; compacting the defaults into
// one wrapper keeps the StorageBackend and SizeTVector typemaps' own messages
// in front of the caller instead of an overload dispatcher's, as for murcko.
%feature("compactdefaultargs") OECluster::take_pairs;
%include "oecluster/Subset.h"

// ============================================================================
// Consensus clustering
// ============================================================================
// All three default their options argument; compacting the defaults into one
// wrapper each keeps the StorageBackend and vector typemaps' own messages in
// front of the caller instead of an overload dispatcher's, as for take_pairs.
// std::vector<int> and std::vector<double> members and returns ride the
// IntVector and DoubleVector templates declared above.
%feature("compactdefaultargs") OECluster::coassociation_distances;
%feature("compactdefaultargs") OECluster::consensus_components;
%feature("compactdefaultargs") OECluster::consensus_strength;

// The allocations below are guarded inside the typemap bodies rather than by
// an %exception: %exception wraps only $action, and an out typemap is emitted
// after that catch ladder has closed (a member getter has none at all), so an
// escaping std::bad_alloc would terminate the interpreter. Both bad_alloc and
// length_error map to SWIG_MemoryError, as in OECLUSTER_GIL_MEMORY_EXCEPTION.

// std_vector.i turns a by-value std::vector return into a tuple, but the
// Python layer is promised an IntVector, so hand back an owned one.
%typemap(out) std::vector<int> consensus_components {
    try {
        $result = SWIG_NewPointerObj(new std::vector<int>($1),
                                     $descriptor(std::vector<int>*),
                                     SWIG_POINTER_OWN);
    } catch (const std::bad_alloc&) {
        SWIG_exception_fail(SWIG_MemoryError, "consensus_components: out of memory");
    } catch (const std::length_error& e) {
        SWIG_exception_fail(SWIG_MemoryError, e.what());
    }
}
// A generated setter assigns a std::vector outside any catch ladder -- an
// out typemap guards only the getter -- so a large assignment would let
// std::bad_alloc escape the extension and terminate the interpreter. These
// two are kernel outputs that nothing assigns to, so drop the setters rather
// than guard them. Named one at a time: ConsensusOptions must stay writable,
// since callers set num_threads and chunk_size on it.
%immutable OECluster::ConsensusStrength::item_consensus;
%immutable OECluster::ConsensusStrength::cluster_consensus;

// The default getter hands back a pointer into the parent struct, which is
// freed as soon as a temporary such as consensus_strength(...).item_consensus
// goes out of scope. %naturalvar would copy but converts to a tuple, so copy
// into an owned DoubleVector instead, scoped to these two members.
%typemap(out) std::vector<double>* item_consensus, std::vector<double>* cluster_consensus {
    try {
        $result = SWIG_NewPointerObj(new std::vector<double>(*$1),
                                     $descriptor(std::vector<double>*),
                                     SWIG_POINTER_OWN);
    } catch (const std::bad_alloc&) {
        SWIG_exception_fail(SWIG_MemoryError, "ConsensusStrength: out of memory");
    } catch (const std::length_error& e) {
        SWIG_exception_fail(SWIG_MemoryError, e.what());
    }
}
%include "oecluster/Consensus.h"

// ============================================================================
// Version macros
// ============================================================================
#define OECLUSTER_VERSION_MAJOR 5
#define OECLUSTER_VERSION_MINOR 26
#define OECLUSTER_VERSION_PATCH 0

// ============================================================================
// Module-level Python convenience code
// ============================================================================
%pythoncode %{
__version__ = "5.26.0"
%}
