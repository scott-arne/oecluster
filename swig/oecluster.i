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
#include "oecluster/PoolMode.h"
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
#include <limits>
#include <optional>
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
// OpenEye's Python bindings use SWIG runtime 4; this module uses runtime 5,
// which oefp shares. SWIG_TypeQuery searches only this module's runtime, so
// it cannot see OpenEye's types, and the conversions below read an object's
// ``this`` themselves. The SWIG descriptor that ``this`` carries decides the
// C++ type, not anything the Python object reports about itself:
// - the descriptor is looked up in the runtime that created it: runtime 5
//   through SWIG's own functions, and runtime 4 through private mirrors of
//   its structures;
// - the pointer is converted through that runtime's cast list;
// - a null descriptor, or one that no readable runtime holds, raises
//   RuntimeError.
// _oecluster_swig_cast, under "Cross-runtime type verification" below, does
// this for each of these conversions and for the grid return path. A
// conversion's isinstance check, where it makes one, runs before
// _oecluster_swig_cast reads the descriptor, and keeps its messages and its
// exception propagation. Overload typechecks never call _oecluster_swig_cast:
// they only choose an overload, and the conversion SWIG then selects calls
// it.
//
// This passes OpenEye and oefp objects between Python and C++ without
// serialization. The macros below generate the boilerplate for each type.

%{
// Returns a new reference to the name an error message gives item's type, or
// NULL with an exception pending for the caller to propagate. The stable ABI
// makes PyTypeObject opaque, so the name is read from the type's __name__
// rather than tp_name.
//
// A metaclass can make __name__ raise or return a non-string. An Exception
// from the read, or a non-string, gives "<unnamed type>": the caller still
// owes its own error, not whatever reading the name raised. Any other
// BaseException, such as KeyboardInterrupt or SystemExit, propagates
// unchanged, because a name lookup must not swallow an interrupt or an exit
// request. RMSD's Python item check applies the same rule.
static PyObject* _oecluster_item_type_name(PyObject* item) {
    PyObject* name = PyObject_GetAttrString((PyObject*)Py_TYPE(item), "__name__");
    if (!name && !PyErr_ExceptionMatches(PyExc_Exception)) {
        return NULL;
    }
    if (!name || !PyUnicode_Check(name)) {
        PyErr_Clear();
        Py_XDECREF(name);
        name = PyUnicode_FromString("<unnamed type>");
    }
    return name;
}

// Sets TypeError "<noun> list item <index> is a <type name><tail>", with
// <index> also stored on the instance as ``_oecluster_item_index``. Call it
// with no exception pending. It returns with one pending, that TypeError or
// whatever stopped it being built, for the caller to propagate.
//
// The attribute marks the error as a conversion's own. cdist renames only a
// marked error, to name the input set (items_named_by_set in
// python/oecluster/_comparisons.py), so an element whose check raised a
// TypeError worded like this one keeps its identity. RMSD's Python item check
// raises a marked copy of this message. A change to the format or the
// attribute must change that module too.
static void _oecluster_list_item_error(const char* noun, Py_ssize_t index, PyObject* item,
                                       const char* tail) {
    PyObject* name = _oecluster_item_type_name(item);
    if (!name) {
        return;
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

// ---- Cross-runtime type verification ----
// _oecluster_swig_cast decides whether an object from another SWIG module
// wraps a given C++ type from the SWIG descriptor its ``this`` carries, not
// from anything the Python object reports about itself. It looks for that
// descriptor by pointer identity in each runtime it can read: this module's
// runtime 5, which oefp shares, and OpenEye's runtime 4. Nothing reached
// through the descriptor is used until a runtime's table is found to hold
// it. The pointer is then converted through that runtime's own cast list,
// which applies whatever adjustment the C++ conversion needs. Runtime 4's
// cast entries are matched by name, as SWIG's own runtime-4 check does.

// Runtime-4 layouts. OpenEye's bindings use SWIG runtime "4", whose
// swig_cast_info is { type, converter, next, prev }; this module compiles
// runtime "5", whose fourth member is an unsigned index instead, so runtime 4
// is read through these private mirrors.
struct _oecluster_swig4_type_info;
struct _oecluster_swig4_cast_info {
    _oecluster_swig4_type_info* type;
    void* (*converter)(void*, int*);
    _oecluster_swig4_cast_info* next;
    _oecluster_swig4_cast_info* prev;
};
struct _oecluster_swig4_type_info {
    const char* name;
    const char* str;
    void* dcast;
    _oecluster_swig4_cast_info* cast;
    void* clientdata;
    int owndata;
};
struct _oecluster_swig4_module_info {
    _oecluster_swig4_type_info** types;
    size_t size;
    _oecluster_swig4_module_info* next;
    void* type_initial;
    void* cast_initial;
    void* clientdata;
};
// The SwigPyObject prefix both runtimes share.
struct _oecluster_swig_object_prefix {
    PyObject_HEAD
    void* ptr;
    void* ty;
};

// Releases the reference _oecluster_swig_cast hands back in *held when the
// scope ends. Like every Py_XDECREF, it must run with the GIL held. Declare
// one, pass &guard.ref, and set ref to null to hand the reference over.
// Copying is deleted because a copy would release the same reference twice.
struct _oecluster_held_ref {
    PyObject* ref = nullptr;
    _oecluster_held_ref() = default;
    _oecluster_held_ref(const _oecluster_held_ref&) = delete;
    _oecluster_held_ref& operator=(const _oecluster_held_ref&) = delete;
    ~_oecluster_held_ref() { Py_XDECREF(ref); }
};

// type.__dict__["__name__"] and type.__dict__["__basicsize__"], read once and
// kept for the process.
static PyObject* _oecluster_type_name_getter = NULL;
static PyObject* _oecluster_type_basicsize_getter = NULL;

// Returns 1 when this_obj's type is really named SwigPyObject and is really
// large enough for _oecluster_swig_object_prefix, 0 when it is not, or -1
// with an exception pending. "Really" means read through type's own
// descriptors, called on the type: attribute lookup on the type goes through
// its metaclass, which can report any name and any size. A Python class
// passes only by being genuinely that large, so reading ptr and ty stays
// within the object.
static int _oecluster_is_swig_object(PyObject* this_obj) {
    if (!_oecluster_type_name_getter || !_oecluster_type_basicsize_getter) {
        PyObject* type_dict = PyObject_GetAttrString((PyObject*)&PyType_Type, "__dict__");
        if (!type_dict) {
            return -1;
        }
        PyObject* name_getter = PyMapping_GetItemString(type_dict, "__name__");
        PyObject* size_getter =
            name_getter ? PyMapping_GetItemString(type_dict, "__basicsize__") : NULL;
        Py_DECREF(type_dict);
        if (!size_getter) {
            Py_XDECREF(name_getter);
            return -1;
        }
        _oecluster_type_name_getter = name_getter;
        _oecluster_type_basicsize_getter = size_getter;
    }
    PyObject* type = (PyObject*)Py_TYPE(this_obj);
    PyObject* name = PyObject_CallMethod(_oecluster_type_name_getter, "__get__", "(O)", type);
    if (!name) {
        return -1;
    }
    bool named = PyUnicode_Check(name) &&
                 PyUnicode_CompareWithASCIIString(name, "SwigPyObject") == 0;
    Py_DECREF(name);
    if (!named) {
        return 0;
    }
    PyObject* size = PyObject_CallMethod(_oecluster_type_basicsize_getter, "__get__", "(O)",
                                         type);
    if (!size) {
        return -1;
    }
    Py_ssize_t basicsize = PyLong_AsSsize_t(size);
    Py_DECREF(size);
    if (basicsize == -1 && PyErr_Occurred()) {
        return -1;
    }
    return basicsize >= (Py_ssize_t)sizeof(_oecluster_swig_object_prefix) ? 1 : 0;
}

// Returns 1 when ty is a descriptor in runtime 5's table, the circular list
// of every runtime-5 module loaded (this one, and oefp's once it is
// imported); otherwise 0.
static int _oecluster_runtime5_holds(void* ty) {
    swig_module_info* head = SWIG_GetModule(0);
    if (!head) {
        return 0;
    }
    swig_module_info* module = head;
    do {
        for (size_t i = 0; i < module->size; ++i) {
            if ((void*)module->types[i] == ty) {
                return 1;
            }
        }
        module = module->next;
    } while (module && module != head);
    return 0;
}

// Runtime 4's module list, once found. It is static memory in an OpenEye
// extension module, which is never unloaded, and runtime-4 modules imported
// later join the same circular list, so it is kept for the process.
static _oecluster_swig4_module_info* _oecluster_swig4_head = NULL;

// Returns 1 with *head set to runtime 4's module list, 0 when runtime 4 is not
// loaded, or -1 with an exception pending. Runtime 4 appears only once an
// OpenEye module is imported, which can be after this module's first
// conversion, so only a success is kept and a failure is retried on the next
// call. swig_runtime_data4 is looked up in sys.modules, never imported: an
// import of an absent module would search the disk for it.
static int _oecluster_runtime4(_oecluster_swig4_module_info** head) {
    if (_oecluster_swig4_head) {
        *head = _oecluster_swig4_head;
        return 1;
    }
    void* pointer = NULL;
    PyObject* runtime = PyMapping_GetItemString(PyImport_GetModuleDict(), "swig_runtime_data4");
    if (runtime) {
        PyObject* capsule = PyObject_GetAttrString(runtime, "type_pointer_capsule");
        Py_DECREF(runtime);
        if (capsule) {
            pointer = PyCapsule_GetPointer(capsule, "swig_runtime_data4.type_pointer_capsule");
            Py_DECREF(capsule);
        }
    }
    if (!pointer) {
        // Each step that returned NULL left an exception. An Exception means
        // runtime 4 is absent or not what SWIG publishes; any other
        // BaseException, such as KeyboardInterrupt, propagates.
        if (!PyErr_ExceptionMatches(PyExc_Exception)) {
            return -1;
        }
        PyErr_Clear();
        return 0;
    }
    _oecluster_swig4_head = (_oecluster_swig4_module_info*)pointer;
    *head = _oecluster_swig4_head;
    return 1;
}

// Returns 1 when ty is a descriptor in one of runtime 4's modules, otherwise 0.
static int _oecluster_runtime4_holds(_oecluster_swig4_module_info* head, void* ty) {
    _oecluster_swig4_module_info* module = head;
    do {
        for (size_t i = 0; i < module->size; ++i) {
            if ((void*)module->types[i] == ty) {
                return 1;
            }
        }
        module = module->next;
    } while (module && module != head);
    return 0;
}

// Returns runtime 4's descriptor with the mangled name name, or NULL when it
// has none. Every module's types table holds the descriptor registered under
// a name, whichever module registered it first, so every match is the same
// descriptor in any import order. Cast entries are not canonicalised that
// way; _oecluster_swig_cast matches them by name.
static _oecluster_swig4_type_info* _oecluster_runtime4_find(_oecluster_swig4_module_info* head,
                                                            const char* name) {
    _oecluster_swig4_module_info* module = head;
    do {
        for (size_t i = 0; i < module->size; ++i) {
            const char* candidate = module->types[i]->name;
            if (candidate && strcmp(candidate, name) == 0) {
                return module->types[i];
            }
        }
        module = module->next;
    } while (module && module != head);
    return NULL;
}

// Decides whether obj wraps target, this module's descriptor for a C++
// pointer type, and if so converts the pointer obj wraps to that type.
// Returns:
//  1: *out holds the converted pointer, and *held a new reference to the
//     verified SwigPyObject that contains it. The caller releases *held only
//     once it no longer uses *out: a ``this`` produced on each read owns its
//     C++ object, and *held may be the only reference keeping it alive.
//  0: obj is not the target type, which includes a ``this`` that is not
//     really a SwigPyObject.
//  2: obj wraps nothing: it has no ``this``, or its ``this`` is a
//     SwigPyObject whose descriptor was found and whose pointer is null.
//  3: obj's ``this`` carries a null descriptor, or one no readable runtime
//     holds.
// -1: an exception is pending for the caller to propagate: reading ``this``
//     raised something other than AttributeError, reading its type failed,
//     looking up runtime 4 raised a BaseException that is not an Exception,
//     or the conversion reported new memory.
// On every return but 1, *out and *held are NULL; 0, 2 and 3 leave no
// exception pending.
//
// No verdict is cached. The same descriptor is asked about different
// targets, such as OEMolBase and then OEMCMolBase, and an answer for one is
// not an answer for the other.
static int _oecluster_swig_cast(PyObject* obj, swig_type_info* target, void** out,
                                PyObject** held) {
    *out = NULL;
    *held = NULL;
    _oecluster_held_ref this_obj;
    this_obj.ref = PyObject_GetAttrString(obj, "this");
    if (!this_obj.ref) {
        if (!PyErr_ExceptionMatches(PyExc_AttributeError)) {
            return -1;
        }
        PyErr_Clear();
        return 2;
    }
    int is_swig_object = _oecluster_is_swig_object(this_obj.ref);
    if (is_swig_object != 1) {
        return is_swig_object;
    }
    const _oecluster_swig_object_prefix* swig_this =
        (const _oecluster_swig_object_prefix*)this_obj.ref;
    void* ty = swig_this->ty;
    void* ptr = swig_this->ptr;
    if (!ty) {
        return 3;
    }
    void* converted = NULL;
    int newmemory = 0;
    if (_oecluster_runtime5_holds(ty)) {
        if (!ptr) {
            return 2;
        }
        swig_cast_info* cast = SWIG_TypeCheckStruct((swig_type_info*)ty, target);
        if (!cast) {
            return 0;
        }
        converted = SWIG_TypeCast(cast, ptr, &newmemory);
    } else {
        _oecluster_swig4_module_info* runtime4 = NULL;
        int loaded = _oecluster_runtime4(&runtime4);
        if (loaded < 0) {
            return -1;
        }
        if (!loaded || !_oecluster_runtime4_holds(runtime4, ty)) {
            return 3;
        }
        if (!ptr) {
            return 2;
        }
        // ty's membership was checked by identity, but a cast entry is matched
        // by its type's name too, because runtime 4 does not canonicalise
        // cast entries. A module that is not the first to register a type adds
        // its casts into that type by name, each still pointing at the
        // module's own copy of the type it casts from. With openeye.oegraphsim
        // imported before openeye.oechem, OEMCMolBase's _p_OEMolWrapper entry
        // is oechem's copy, not the descriptor an OEMol's ``this`` carries.
        // SWIG's own runtime-4 check compares names for the same reason.
        // Matching by name admits no type the list does not name: OEMCMolBase's
        // list has no _p_OEGraphMolWrapper entry under any name, in either
        // import order. ty is a genuine descriptor once found, and each
        // entry's type is runtime 4's own, so both names are safe to read.
        //
        // Runtime 4's own type checks move the cast they find to the front of
        // this list. They and this walk both run with the GIL held, so the
        // walk never sees the list mid-change.
        _oecluster_swig4_type_info* target4 = _oecluster_runtime4_find(runtime4, target->name);
        const char* ty_name = ((const _oecluster_swig4_type_info*)ty)->name;
        _oecluster_swig4_cast_info* cast = target4 ? target4->cast : NULL;
        for (; cast; cast = cast->next) {
            const char* cast_name = cast->type->name;
            if ((void*)cast->type == ty ||
                (ty_name && cast_name && strcmp(cast_name, ty_name) == 0)) {
                break;
            }
        }
        if (!cast) {
            return 0;
        }
        converted = cast->converter ? cast->converter(ptr, &newmemory) : ptr;
    }
    if (newmemory) {
        // The new object is left unreleased: this code knows neither its type
        // nor how to free it.
        PyObject* name = _oecluster_item_type_name(obj);
        if (name) {
            PyErr_Format(PyExc_RuntimeError,
                         "converting %U reported new memory, which oecluster cannot release",
                         name);
            Py_DECREF(name);
        }
        return -1;
    }
    *out = converted;
    *held = this_obj.ref;
    this_obj.ref = nullptr;
    return 1;
}

// Returns a new reference to the names of the swig_runtime_data<N> modules
// in sys.modules, sorted and joined with ", ", or "none" when there are none.
// Returns NULL with an exception pending if the list cannot be built. After
// any failure only releases run before it returns, so no other call runs
// with that exception pending.
static PyObject* _oecluster_loaded_runtimes() {
    PyObject* keys = PyMapping_Keys(PyImport_GetModuleDict());
    if (!keys) {
        return NULL;
    }
    PyObject* prefix = PyUnicode_FromString("swig_runtime_data");
    if (!prefix) {
        Py_DECREF(keys);
        return NULL;
    }
    PyObject* runtimes = PyList_New(0);
    if (!runtimes) {
        Py_DECREF(prefix);
        Py_DECREF(keys);
        return NULL;
    }
    Py_ssize_t count = PyList_Size(keys);
    for (Py_ssize_t i = 0; i < count; ++i) {
        PyObject* key = PyList_GetItem(keys, i);
        Py_ssize_t match = PyUnicode_Check(key)
                               ? PyUnicode_Tailmatch(key, prefix, 0, PY_SSIZE_T_MAX, -1)
                               : 0;
        if (match < 0 || (match == 1 && PyList_Append(runtimes, key) < 0)) {
            Py_DECREF(runtimes);
            Py_DECREF(prefix);
            Py_DECREF(keys);
            return NULL;
        }
    }
    Py_DECREF(prefix);
    Py_DECREF(keys);
    if (PyList_Sort(runtimes) < 0) {
        Py_DECREF(runtimes);
        return NULL;
    }
    if (PyList_Size(runtimes) == 0) {
        Py_DECREF(runtimes);
        return PyUnicode_FromString("none");
    }
    PyObject* separator = PyUnicode_FromString(", ");
    if (!separator) {
        Py_DECREF(runtimes);
        return NULL;
    }
    PyObject* result = PyUnicode_Join(separator, runtimes);
    Py_DECREF(separator);
    Py_DECREF(runtimes);
    return result;
}

// Sets RuntimeError for an item _oecluster_swig_cast answered 3 for, as
// "<noun> list item is a <type name> whose SWIG type is ..." from a list
// conversion, or "argument is a <type name> whose SWIG type is ..." when noun
// is NULL. Call it with no exception pending. It returns with one pending,
// that RuntimeError or whatever stopped it being built, for the caller to
// propagate. The message carries no index: cdist renames only a marked
// TypeError, so an index would reach the caller as one into
// items_a + items_b.
static void _oecluster_unreadable_runtime_error(const char* noun, PyObject* item) {
    PyObject* name = _oecluster_item_type_name(item);
    if (!name) {
        return;
    }
    PyObject* loaded = _oecluster_loaded_runtimes();
    if (!loaded) {
        Py_DECREF(name);
        return;
    }
    if (noun) {
        PyErr_Format(PyExc_RuntimeError,
                     "%s list item is a %U whose SWIG type is in no SWIG runtime oecluster "
                     "can read (loaded: %U)",
                     noun, name, loaded);
    } else {
        PyErr_Format(PyExc_RuntimeError,
                     "argument is a %U whose SWIG type is in no SWIG runtime oecluster can "
                     "read (loaded: %U)",
                     name, loaded);
    }
    Py_DECREF(loaded);
    Py_DECREF(name);
}

// The fallback of the OE_CROSS_RUNTIME_*_TYPEMAPS macros, for an argument
// SWIG_ConvertPtr has not converted with this module's descriptor. It returns
// _oecluster_swig_cast's code, with *argp and *held set as that describes.
// It lives here, not in the macros, so that it is compiled in every build:
// a macro's body is compiled only where a wrapped function expands it.
static int _oecluster_foreign_ref(PyObject* input, swig_type_info* descriptor, void** argp,
                                  PyObject** held) {
    return _oecluster_swig_cast(input, descriptor, argp, held);
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

// Returns _oecluster_swig_cast's code for obj against target, this module's
// descriptor for OEFP::OEFPBatch *, which the typemap passes in. A high-level
// OEFPBatch wraps the native batch in ``_native``, so the native batch is the
// object checked; any other object is checked itself. ``_native`` is read
// under the rule for ``this``: AttributeError means the object wraps nothing,
// 2, and any other exception is -1, left pending.
//
// On 1, *held is the helper's own reference to the verified ``this``. A
// ``this`` that owns its C++ batch keeps it alive by itself, so ``_native``
// is released here, whatever the code.
static int _oecluster_extract_oefp_batch_ptr(PyObject* obj, swig_type_info* target, void** out,
                                             PyObject** held) {
    *out = NULL;
    *held = NULL;
    if (_oecluster_load_oefp_batch_types() < 0) {
        return -1;
    }
    int is_high_level = PyObject_IsInstance(obj, _oecluster_oefp_batch_type);
    if (is_high_level < 0) {
        return -1;
    }
    if (!is_high_level) {
        return _oecluster_swig_cast(obj, target, out, held);
    }
    PyObject* native_obj = PyObject_GetAttrString(obj, "_native");
    if (!native_obj) {
        if (!PyErr_ExceptionMatches(PyExc_AttributeError)) {
            return -1;
        }
        PyErr_Clear();
        return 2;
    }
    int code = _oecluster_swig_cast(native_obj, target, out, held);
    Py_DECREF(native_obj);
    return code;
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
// exception is then lost, and the in typemap goes by the second answer: it
// refuses the object on a no, and on a yes passes it to the descriptor check,
// which decides the type (see the cross-runtime comment above). Only a
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
    // Nothing is deleted or written through until the new object's ``this``
    // is verified as an OEScalarGrid. The write then goes to the exact
    // SwigPyObject that was verified, not to whatever a second read of
    // ``this`` returns. A newly constructed object owns its pointer, so once
    // grid is stored there, that object frees it. A check that raised keeps
    // its own exception, as at every other conversion.
    _oecluster_held_ref verified_this;
    void* old_ptr = NULL;
    int cast = _oecluster_swig_cast(oe_grid, SWIGTYPE_p_OESystem__OEScalarGrid, &old_ptr,
                                    &verified_this.ref);
    if (cast != 1) {
        delete grid;
        // Nothing the new object holds has been touched, so releasing it
        // frees whatever it owns.
        Py_DECREF(oe_grid);
        if (cast != -1) {
            PyErr_SetString(PyExc_RuntimeError,
                            "oecluster could not verify the OEScalarGrid it created");
        }
        return NULL;
    }
    delete static_cast<OESystem::OEScalarGrid*>(old_ptr);
    ((_oecluster_swig_object_prefix*)verified_this.ref)->ptr = grid;
    return oe_grid;
}
%}

// Descriptors for OpenEye pointer types, for use as _oecluster_swig_cast's
// targets. Declaring them puts them in this module's type table whether or
// not a wrapped function names them. _swig_cast_verdict looks them up by
// name, and _oecluster_swig_cast searches runtime 4 by the same names.
%types(OEChem::OEMolBase *, OEChem::OEMCMolBase *, OEBio::OEDesignUnit *,
       OESystem::OEScalarGrid *);

%rename(_swig_cast_verdict) oecluster_swig_cast_verdict;
%feature("docstring") oecluster_swig_cast_verdict
"For tests only. Returns the code oecluster's cross-runtime type check gives obj
against the descriptor named target_name: 1 converts, 0 does not, 2 wraps
nothing, 3 is in no readable SWIG runtime.";
%inline %{
// For tests only: the helper's verdict for obj against the descriptor named
// target_name, with the converted pointer and the held reference discarded.
PyObject* oecluster_swig_cast_verdict(PyObject* obj, const char* target_name) {
    swig_type_info* target = SWIG_MangledTypeQuery(target_name);
    if (!target) {
        PyErr_Format(PyExc_ValueError, "no descriptor named %s in oecluster", target_name);
        return NULL;
    }
    void* out = NULL;
    PyObject* held = NULL;
    int code = _oecluster_swig_cast(obj, target, &out, &held);
    Py_XDECREF(held);
    if (code < 0) {
        return NULL;
    }
    return PyLong_FromLong(code);
}
%}

// ============================================================================
// Typemap generator macros
// ============================================================================

// Generate const-ref and non-const-ref typemaps for a cross-runtime OpenEye type.
// Each in typemap tries SWIG_ConvertPtr with this module's descriptor first.
// Failing that, an object the TAG check accepts goes to
// _oecluster_foreign_ref, and its descriptor decides the cast: an object that
// is not the type, or wraps nothing, gets ERR_MSG, and one whose runtime
// cannot be read gets the unreadable-runtime RuntimeError. The verified
// ``this`` is kept until the matching freearg releases it, once the call
// has returned or failed.
// The nullable-pointer macro below follows the same rule.
// CPP_TYPE: fully qualified C++ type (e.g., OEChem::OEMolBase)
// TAG:      DEFINE_OE_TYPE_CHECKER suffix (e.g., oemolbase); the in typemaps
//           call _oecluster_isinstance_TAG and the typechecks
//           _oecluster_may_be_TAG
// ERR_MSG:  error message on type mismatch
%define OE_CROSS_RUNTIME_REF_TYPEMAPS(CPP_TYPE, TAG, ERR_MSG)

%typemap(in) const CPP_TYPE& (void *argp = 0, int res = 0, PyObject* held = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        int is_instance = _oecluster_isinstance_##TAG($input);
        if (is_instance < 0) {
            SWIG_fail;
        }
        if (is_instance) {
            int cast = _oecluster_foreign_ref($input, $descriptor, &argp, &held);
            if (cast == 1) {
                res = SWIG_OK;
            } else if (cast == 3) {
                _oecluster_unreadable_runtime_error(NULL, $input);
                SWIG_fail;
            } else if (cast < 0) {
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

%typemap(freearg) const CPP_TYPE& {
    Py_XDECREF(held$argnum);
}

%typemap(typecheck, precedence=10) const CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : _oecluster_may_be_##TAG($input) ? 1 : 0;
}

%typemap(in) CPP_TYPE& (void *argp = 0, int res = 0, PyObject* held = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        int is_instance = _oecluster_isinstance_##TAG($input);
        if (is_instance < 0) {
            SWIG_fail;
        }
        if (is_instance) {
            int cast = _oecluster_foreign_ref($input, $descriptor, &argp, &held);
            if (cast == 1) {
                res = SWIG_OK;
            } else if (cast == 3) {
                _oecluster_unreadable_runtime_error(NULL, $input);
                SWIG_fail;
            } else if (cast < 0) {
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

%typemap(freearg) CPP_TYPE& {
    Py_XDECREF(held$argnum);
}

%typemap(typecheck, precedence=10) CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : _oecluster_may_be_##TAG($input) ? 1 : 0;
}

%enddef

// Generate nullable-pointer typemaps (accepts None) for a cross-runtime type.
%define OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(CPP_TYPE, TAG, ERR_MSG)

%typemap(in) const CPP_TYPE* (void *argp = 0, int res = 0, PyObject* held = 0) {
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
                int cast = _oecluster_foreign_ref($input, $descriptor, &argp, &held);
                if (cast == 1) {
                    res = SWIG_OK;
                } else if (cast == 3) {
                    _oecluster_unreadable_runtime_error(NULL, $input);
                    SWIG_fail;
                } else if (cast < 0) {
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

%typemap(freearg) const CPP_TYPE* {
    Py_XDECREF(held$argnum);
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
%typemap(in) const std::vector<OEChem::OEMolBase*>&
        (std::vector<OEChem::OEMolBase*> temp, PyObject* held = 0, PyObject* this_refs = 0) {
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
    //
    // this_refs holds, for the same span, the ``this`` each pointer was
    // verified and converted from. Holding the element is not enough: a
    // ``this`` produced on each read owns its C++ object, and nothing but
    // the reference _oecluster_swig_cast hands back keeps that object alive.
    held = PyList_AsTuple($input);
    if (!held) {
        SWIG_fail;
    }
    Py_ssize_t size = PyTuple_GET_SIZE(held);
    this_refs = PyTuple_New(size);
    if (!this_refs) {
        SWIG_fail;
    }
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
        // An overridden __class__ can satisfy the check above, so the
        // descriptor the item's ``this`` carries decides the cast. An
        // OEMolBase instance whose ``this`` is missing or null, such as one
        // made by ``__new__`` without ``__init__``, passes the check above,
        // so it gets its own message rather than being called a
        // non-OEMolBase.
        _oecluster_held_ref verified_this;
        void* ptr = NULL;
        int cast = _oecluster_swig_cast(item, $descriptor(OEChem::OEMolBase *), &ptr,
                                        &verified_this.ref);
        if (cast != 1) {
            if (cast == 0) {
                _oecluster_molecule_item_error(i, item, ", not an OEMolBase");
            } else if (cast == 2) {
                _oecluster_molecule_item_error(i, item, " that wraps no molecule");
            } else if (cast == 3) {
                _oecluster_unreadable_runtime_error("molecule", item);
            }
            SWIG_fail;
        }
        // PyTuple_SetItem steals the reference even when it fails, so the
        // guard gives it up first.
        PyObject* verified = verified_this.ref;
        verified_this.ref = nullptr;
        if (PyTuple_SetItem(this_refs, i, verified) < 0) {
            SWIG_fail;
        }
        // The reserve above means this cannot reallocate today; the catch
        // below keeps that from being load-bearing.
        try {
            temp.push_back(static_cast<OEChem::OEMolBase*>(ptr));
        } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a molecule list")
    }
    $1 = &temp;
}

// Py_XDECREF needs the GIL. This relies on every %exception handler having
// reacquired it before the wrapper reaches freearg, on the fail: path as well
// as on success.
%typemap(freearg) const std::vector<OEChem::OEMolBase*>& {
    Py_XDECREF(held$argnum);
    Py_XDECREF(this_refs$argnum);
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
        // The cast target is OEMCMolBase, not OEMol. OpenEye's Python bindings
        // hand out an OEMolWrapper*, which derives from OEMCMolBase but is not
        // an OEChem::OEMol -- that type is a bare handle holding an
        // OEOwnedPtr<OEMCMolBase>. Casting to OEMol and copy-constructing reads
        // the wrapper's first member as if it were that owned pointer and
        // segfaults. OEMol(const OEMCMolBase&) is the documented converting
        // constructor for exactly this, and it preserves conformers.
        //
        // Whether the item wraps an OEMCMolBase is decided by the descriptor
        // its ``this`` carries, not by the check above, which an overridden
        // __class__ can satisfy: an object whose __class__ answers OEMol but
        // which wraps an OEGraphMol is refused here. verified_this keeps the
        // object ptr points to alive until the copy is made.
        _oecluster_held_ref verified_this;
        void* ptr = NULL;
        int cast = _oecluster_swig_cast(item, $descriptor(OEChem::OEMCMolBase *), &ptr,
                                        &verified_this.ref);
        if (cast != 1) {
            if (cast == 0) {
                _oecluster_molecule_item_error(i, item, ", not an OEMol");
            } else if (cast == 2) {
                _oecluster_molecule_item_error(i, item, " that wraps no molecule");
            } else if (cast == 3) {
                _oecluster_unreadable_runtime_error("molecule", item);
            }
            SWIG_fail;
        }
        // Make a COPY wrapped in shared_ptr (we do not own the Python object).
        auto* mcmol = static_cast<OEChem::OEMCMolBase*>(ptr);
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
    // hole an emptied list opens. Such a change would also have to keep each
    // verified ``this`` for the call, as that conversion does: here each one
    // is released once its element is copied.
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
        // Kept so that an exception it raises still propagates. Its answer no
        // longer chooses the copy: an overridden __class__ can make it true
        // for an object that wraps an OEGraphMol.
        int is_multi = _oecluster_isinstance_oemcmolbase(item);
        if (is_multi < 0) {
            SWIG_fail;
        }
        _oecluster_held_ref verified_this;
        void* ptr = NULL;
        int cast = _oecluster_swig_cast(item, $descriptor(OEChem::OEMolBase *), &ptr,
                                        &verified_this.ref);
        if (cast != 1) {
            if (cast == 0) {
                _oecluster_molecule_item_error(i, item, ", not an OEMolBase");
            } else if (cast == 2) {
                _oecluster_molecule_item_error(i, item, " that wraps no molecule");
            } else if (cast == 3) {
                _oecluster_unreadable_runtime_error("molecule", item);
            }
            SWIG_fail;
        }
        // The descriptor chooses the copy. Reading ``this`` again can give a
        // different object, so each pointer is copied only while the
        // reference verified alongside it is kept.
        _oecluster_held_ref verified_multi_this;
        void* multi_ptr = NULL;
        int multi = _oecluster_swig_cast(item, $descriptor(OEChem::OEMCMolBase *), &multi_ptr,
                                         &verified_multi_this.ref);
        if (multi == 2) {
            _oecluster_molecule_item_error(i, item, " that wraps no molecule");
            SWIG_fail;
        }
        if (multi == 3) {
            _oecluster_unreadable_runtime_error("molecule", item);
            SWIG_fail;
        }
        if (multi < 0) {
            SWIG_fail;
        }
        // The copy runs before $action, outside every %exception handler, so
        // it is guarded here by OECLUSTER_TYPEMAP_ALLOCATION_CATCH.
        try {
            if (multi == 1) {
                temp.push_back(std::make_shared<OEChem::OEMol>(
                    *static_cast<OEChem::OEMCMolBase*>(multi_ptr)));
            } else {
                temp.push_back(std::make_shared<OEChem::OEMol>(
                    *static_cast<OEChem::OEMolBase*>(ptr)));
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
        // As in the strict OEMol list conversion, the descriptor decides the
        // cast, and verified_this keeps the object alive until it is copied.
        _oecluster_held_ref verified_this;
        void* ptr = NULL;
        int cast = _oecluster_swig_cast(item, $descriptor(OEBio::OEDesignUnit *), &ptr,
                                        &verified_this.ref);
        if (cast != 1) {
            if (cast == 0) {
                _oecluster_list_item_error("design unit", i, item, ", not an OEDesignUnit");
            } else if (cast == 2) {
                _oecluster_list_item_error("design unit", i, item, " that wraps no design unit");
            } else if (cast == 3) {
                _oecluster_unreadable_runtime_error("design unit", item);
            }
            SWIG_fail;
        }
        // Make a COPY wrapped in shared_ptr (we do not own the Python object)
        auto* du = static_cast<OEBio::OEDesignUnit*>(ptr);
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

// ---- std::function<void(size_t, size_t)> progress callback helpers ----
// The progress typemaps themselves sit with the PDistOptions redeclaration.
%{
// The callable a progress callback typemap wraps. Named, rather than a lambda,
// so the progress getter can recognise it with std::function::target and
// hand back the Python callable it holds.
struct _oecluster_progress_callback {
    std::shared_ptr<PyObject> callable;
    void operator()(size_t completed, size_t total) const {
        PyGILState_STATE gstate = PyGILState_Ensure();
        PyObject* result = PyObject_CallFunction(callable.get(), "nn",
            (Py_ssize_t)completed, (Py_ssize_t)total);
        Py_XDECREF(result);
        if (PyErr_Occurred()) PyErr_Clear();
        PyGILState_Release(gstate);
    }
};

// Builds the std::function for a callable, under the caller's allocation
// guard; the shared_ptr's deleter drops the reference with the GIL held.
static std::function<void(size_t, size_t)> _oecluster_make_progress(PyObject* callable) {
    Py_INCREF(callable);
    std::shared_ptr<PyObject> guard(callable, [](PyObject* p) {
        PyGILState_STATE gstate = PyGILState_Ensure();
        Py_DECREF(p);
        PyGILState_Release(gstate);
    });
    return _oecluster_progress_callback{guard};
}
%}

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

%{
// The memory-pool mode the application set through OpenEye's Python layer,
// read while the wrapper still holds the GIL. Looked up in sys.modules rather
// than imported: importing oecluster does not import openeye.oechem, and an
// application that never loaded it cannot have set a mode through it, so the
// native mode stands. A failed read answers 0, an unsafe mode: a call wrongly
// kept serial is slow, while one wrongly run in parallel crashes. A result
// that is not a plain int within unsigned int's range also answers 0. Only a
// cleanly missing module answers "no value".
static std::optional<unsigned int> _oecluster_caller_pool_mode() {
    PyObject* modules = PyImport_GetModuleDict();  // borrowed
    if (modules == nullptr) {
        PyErr_Clear();
        return 0u;
    }
    PyObject* key = PyUnicode_FromString("openeye.oechem");
    if (key == nullptr) {
        PyErr_Clear();
        return 0u;
    }
    // PyDict_GetItemWithError tells a missing key (no error set) from a
    // failed lookup (error set); only the first means "never loaded".
    PyObject* oechem = PyDict_GetItemWithError(modules, key);  // borrowed
    Py_DECREF(key);
    if (oechem == nullptr) {
        if (PyErr_Occurred()) {
            PyErr_Clear();
            return 0u;
        }
        return std::nullopt;
    }
    PyObject* value = PyObject_CallMethod(oechem, "OEGetMemPoolMode", nullptr);
    if (value == nullptr) {
        PyErr_Clear();
        return 0u;
    }
    // Exact, because bool is an int subclass and True would read as 1, which
    // is Mutexed and thread-safe.
    if (!PyLong_CheckExact(value)) {
        Py_DECREF(value);
        return 0u;
    }
    const unsigned long mode = PyLong_AsUnsignedLong(value);
    Py_DECREF(value);
    if (PyErr_Occurred()) {
        PyErr_Clear();
        return 0u;
    }
    // unsigned long is wider than unsigned int on LP64, where the cast below
    // would turn 2**32 + 1 into 1.
    if (mode > std::numeric_limits<unsigned int>::max()) {
        return 0u;
    }
    return static_cast<unsigned int>(mode);
}
%}

%rename(_caller_pool_mode) oecluster_caller_pool_mode;
%feature("docstring") oecluster_caller_pool_mode
"For tests only. Returns the memory-pool mode oecluster reads from
openeye.oechem.OEGetMemPoolMode(), 0 when that read fails or gives anything but
a plain int (not a bool or other int subclass) within unsigned int's range, or
None when openeye.oechem was never imported.";
%inline %{
// For tests only: the reader's answer, as a Python int or None.
PyObject* oecluster_caller_pool_mode() {
    const std::optional<unsigned int> mode = _oecluster_caller_pool_mode();
    if (!mode) {
        Py_RETURN_NONE;
    }
    return PyLong_FromUnsignedLong(*mode);
}
%}

// ============================================================================
// GIL release for the long-running entry points (overrides the general
// handler above)
//
// Every override releases the GIL around $action, unless the memory-pool mode
// is unsafe (below), and restores it on each exit path. Most of them then
// share one catch ladder, differing only in the name in the fallback message;
// spelling that out by hand for each one invited the ladders to drift apart,
// and one of them already had. `NAME` is SWIG preprocessor stringification,
// which is what keeps each message naming its own function. The entry points
// that can plausibly exhaust memory use OECLUSTER_GIL_MEMORY_EXCEPTION
// instead, which adds a pair of catch clauses.
//
// Under a resolved memory-pool mode that is not thread-safe, the call keeps
// the GIL. OpenEye then requires single-threaded use, and holding the GIL
// serialises the call against every other Python thread. Another thread runs
// only while a progress callback runs Python code, and the call's OpenEye
// work, then all on this thread, waits for it. The GIL cannot hold back an
// OpenEye call another thread is running, because OpenEye's wrappers release
// the GIL around their work (measured on 2026-10-24); OpenEye's own rule under
// these modes forbids such a call. The scope records the caller's mode for
// the native worker gate, and tells it that a progress callback can run only
// on this thread. The inner braces keep the initialized locals out of the
// scope that SWIG's earlier `goto fail` jumps cross. The scope is destroyed on
// every exit, including SWIG_exception's jump to `fail`.
// ============================================================================
%define OECLUSTER_GIL_EXCEPTION(QUALIFIED, NAME)
%exception QUALIFIED {
    {
        const std::optional<unsigned int> _oecluster_python_mode =
            _oecluster_caller_pool_mode();
        const bool _oecluster_keep_gil = !OECluster::pool_is_thread_safe(
            OECluster::resolve_pool_mode(_oecluster_python_mode));
        OECluster::CallerPoolModeScope _oecluster_pool_scope(
            _oecluster_python_mode, _oecluster_keep_gil);
        PyThreadState* _oecluster_saved =
            _oecluster_keep_gil ? nullptr : PyEval_SaveThread();
        try {
            $action
        } catch (const OECluster::OEClusterError& e) {
            if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
            SWIG_exception(SWIG_RuntimeError, e.what());
        } catch (const std::exception& e) {
            if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
            SWIG_exception(SWIG_RuntimeError, e.what());
        } catch (...) {
            if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
            SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception in " `NAME`);
        }
        if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
    }
}
%enddef

// The same ladder, except that std::bad_alloc and std::length_error (a
// container asked to exceed its max_size()) surface as SWIG_MemoryError rather
// than the blanket SWIG_RuntimeError, so a caller can tell an exhausted
// machine from a refused input.
%define OECLUSTER_GIL_MEMORY_EXCEPTION(QUALIFIED, NAME)
%exception QUALIFIED {
    {
        const std::optional<unsigned int> _oecluster_python_mode =
            _oecluster_caller_pool_mode();
        const bool _oecluster_keep_gil = !OECluster::pool_is_thread_safe(
            OECluster::resolve_pool_mode(_oecluster_python_mode));
        OECluster::CallerPoolModeScope _oecluster_pool_scope(
            _oecluster_python_mode, _oecluster_keep_gil);
        PyThreadState* _oecluster_saved =
            _oecluster_keep_gil ? nullptr : PyEval_SaveThread();
        try {
            $action
        } catch (const OECluster::OEClusterError& e) {
            if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
            SWIG_exception(SWIG_RuntimeError, e.what());
        } catch (const std::bad_alloc&) {
            if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
            SWIG_exception(SWIG_MemoryError, `NAME` ": out of memory");
        } catch (const std::length_error& e) {
            if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
            SWIG_exception(SWIG_MemoryError, e.what());
        } catch (const std::exception& e) {
            if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
            SWIG_exception(SWIG_RuntimeError, e.what());
        } catch (...) {
            if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
            SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception in " `NAME`);
        }
        if (_oecluster_saved) PyEval_RestoreThread(_oecluster_saved);
    }
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

// For the descriptor entry points, which compute descriptors through oefp
// threads and run under the general ladder with the GIL held: the scope lets
// calculate_descriptor_batch see the caller's mode in a statically linked
// build. Not added to the general %exception, which would put a sys.modules
// read in front of every wrapped call.
%define OECLUSTER_POOL_SCOPED_EXCEPTION(QUALIFIED)
%exception QUALIFIED {
    {
        OECluster::CallerPoolModeScope _oecluster_pool_scope(
            _oecluster_caller_pool_mode());
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
}
%enddef

OECLUSTER_POOL_SCOPED_EXCEPTION(OECluster::DescriptorComparison::DescriptorComparison)
OECLUSTER_POOL_SCOPED_EXCEPTION(OECluster::descriptor_statistics)
OECLUSTER_POOL_SCOPED_EXCEPTION(OECluster::descriptor_excluded_indices)

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

// Internal capability read by the worker gate (oecluster/PoolMode.h), not
// Python API.
%ignore RunsOpenEyeCode;

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

%typemap(in) const OEFP::OEFPBatch& (void *argp = 0, PyObject* held = 0) {
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
    // An overridden __class__ can satisfy the check above, so the descriptor
    // decides the cast. The verified ``this`` is kept until the freearg
    // below, because the call reads the batch through the pointer.
    int cast = _oecluster_extract_oefp_batch_ptr($input, $descriptor(OEFP::OEFPBatch *), &argp,
                                                 &held);
    if (cast == 0) {
        SWIG_exception_fail(SWIG_TypeError, "Expected an oefp.OEFPBatch.");
    }
    if (cast == 2) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null OEFPBatch reference.");
    }
    if (cast == 3) {
        _oecluster_unreadable_runtime_error(NULL, $input);
        SWIG_fail;
    }
    if (cast != 1) {
        SWIG_fail;
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

// Same GIL precondition as the OEMolBase list's freearg.
%typemap(freearg) const OEFP::OEFPBatch& {
    Py_XDECREF(held$argnum);
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
// SWIG 4.5.1 applies no %exception to member-variable accessors, so a
// std::bad_alloc from copying a std::function in the progress setter would
// escape the extension and terminate the interpreter; every copy below sits
// in its own typemap's catch ladder instead. %naturalvar passes the member by
// const reference, so the getter copies nothing: the out typemap reads the
// callable through the reference rather than returning a heap copy SWIG has
// no destructor for, and the setter's const-reference in typemap applies.
// The setter's one copy is the assignment in memberin, where it can be
// guarded. %naturalvar acts only on types SWIG knows as classes, and nothing
// else in the interface declares std::function, hence the declaration below.
// It is parse-only, so it never reaches the wrapper, and the empty-name
// %template() registers the type without generating a proxy class. The
// typemaps and %naturalvar are named for the member ``progress`` so that they
// claim no other parameter or member of this type, since their messages speak
// of progress callbacks.
namespace std {
  template <typename> class function;
  template <typename RET, typename... ARGS> class function<RET(ARGS...)> {};
}
%template() std::function<void(size_t, size_t)>;

%naturalvar OECluster::PDistOptions::progress;
%naturalvar OECluster::CDistOptions::progress;

%typemap(in) const std::function<void(size_t, size_t)>& progress
        (std::function<void(size_t, size_t)> temp) {
    if ($input == Py_None) {
        temp = nullptr;
    } else if (PyCallable_Check($input)) {
        // _oecluster_make_progress takes a reference, then makes two
        // allocations: the shared_ptr its control block, and the
        // std::function its copy of the functor when that is stored on the
        // heap. Neither leaks the reference if it throws. shared_ptr calls
        // the deleter itself when the control block cannot be allocated, and
        // a throwing std::function construction unwinds through guard's
        // destructor, before the catch runs.
        try {
            temp = _oecluster_make_progress($input);
        } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("converting a progress callback")
    } else {
        SWIG_exception_fail(SWIG_TypeError, "Expected callable or None for progress callback");
    }
    $1 = &temp;
}

// A memberin typemap replaces SWIG's whole assignment, including its
// ``if (arg1)``, and SWIG_ConvertPtr accepts None as a null self, so the
// check is restored here.
%typemap(memberin) std::function<void(size_t, size_t)> progress {
    if ($self) {
        try {
            $1 = *$input;
        } OECLUSTER_TYPEMAP_ALLOCATION_CATCH("setting a progress callback")
    }
}

%typemap(out) const std::function<void(size_t, size_t)>& progress {
    if (!*$1) {
        $result = Py_None;
        Py_INCREF($result);
    } else if (const _oecluster_progress_callback* held =
                   $1->target<_oecluster_progress_callback>()) {
        $result = held->callable.get();
        Py_INCREF($result);
    } else {
        SWIG_exception_fail(SWIG_RuntimeError,
                            "progress holds a C++ callback that has no Python form");
    }
}

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
#define OECLUSTER_VERSION_PATCH 1

// ============================================================================
// Module-level Python convenience code
// ============================================================================
%pythoncode %{
__version__ = "5.26.1"
%}
