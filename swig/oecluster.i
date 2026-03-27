/* -*- mode: c++ -*- */
%module _oecluster

/* --------------------------------------------------------------------------
 * 1. C++ header includes (compiled into the wrapper, not parsed by SWIG)
 * -------------------------------------------------------------------------- */
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

using namespace OECluster;
%}

/* --------------------------------------------------------------------------
 * 2. Cross-module SWIG type resolution helpers
 * -------------------------------------------------------------------------- */
%{
struct _SwigPyObjectCompat {
    PyObject_HEAD
    void *ptr;
};

static PyObject* _oecluster_oe_molbase_type = NULL;
static PyObject* _oecluster_oe_mol_type = NULL;
static PyObject* _oecluster_oe_du_type = NULL;

static bool _oecluster_is_oemolbase(PyObject* obj) {
    if (!_oecluster_oe_molbase_type) {
        PyObject* mod = PyImport_ImportModule("openeye.oechem");
        if (mod) {
            _oecluster_oe_molbase_type = PyObject_GetAttrString(mod, "OEMolBase");
            Py_DECREF(mod);
        }
        if (!_oecluster_oe_molbase_type) return false;
    }
    return PyObject_IsInstance(obj, _oecluster_oe_molbase_type) == 1;
}

static bool _oecluster_is_oemol(PyObject* obj) {
    if (!_oecluster_oe_mol_type) {
        PyObject* mod = PyImport_ImportModule("openeye.oechem");
        if (mod) {
            _oecluster_oe_mol_type = PyObject_GetAttrString(mod, "OEMol");
            Py_DECREF(mod);
        }
        if (!_oecluster_oe_mol_type) return false;
    }
    return PyObject_IsInstance(obj, _oecluster_oe_mol_type) == 1;
}

static bool _oecluster_is_oedesignunit(PyObject* obj) {
    if (!_oecluster_oe_du_type) {
        PyObject* mod = PyImport_ImportModule("openeye.oebio");
        if (mod) {
            _oecluster_oe_du_type = PyObject_GetAttrString(mod, "OEDesignUnit");
            Py_DECREF(mod);
        }
        if (!_oecluster_oe_du_type) return false;
    }
    return PyObject_IsInstance(obj, _oecluster_oe_du_type) == 1;
}

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
%}

/* --------------------------------------------------------------------------
 * 3. Typemap: const std::vector<OEChem::OEMolBase*>&
 *    Used by FingerprintMetric and SuperposeMetric (mol overload)
 * -------------------------------------------------------------------------- */
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

/* --------------------------------------------------------------------------
 * 4. Typemap: const std::vector<std::shared_ptr<OEChem::OEMol>>&
 *    Used by ROCSMetric
 * -------------------------------------------------------------------------- */
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
        /* Make a COPY wrapped in shared_ptr (we do not own the Python object) */
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

/* --------------------------------------------------------------------------
 * 5. Typemap: const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>&
 *    Used by SuperposeMetric (DU overload)
 * -------------------------------------------------------------------------- */
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
        /* Make a COPY wrapped in shared_ptr (we do not own the Python object) */
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

/* --------------------------------------------------------------------------
 * 6. Typemap: std::function<void(size_t, size_t)> progress callback
 *    Used by PDistOptions::progress
 * -------------------------------------------------------------------------- */
%typemap(in) std::function<void(size_t, size_t)> {
    if ($input == Py_None) {
        $1 = nullptr;
    } else if (PyCallable_Check($input)) {
        PyObject* callback = $input;
        Py_INCREF(callback);
        /* Use shared_ptr guard to ensure Py_DECREF when lambda is destroyed */
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
    /* Prevent copy-on-scope-exit; the shared_ptr guard handles Py_DECREF */
}

/* --------------------------------------------------------------------------
 * 7. STL support and common SWIG library includes
 * -------------------------------------------------------------------------- */
%include "std_string.i"
%include "std_vector.i"
%include "stdint.i"
%include "exception.i"

/* Instantiate vector templates used by the Python layer */
%template(StringVector) std::vector<std::string>;

/* --------------------------------------------------------------------------
 * 8. Exception handling -- general
 * -------------------------------------------------------------------------- */
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

/* --------------------------------------------------------------------------
 * 9. GIL release for pdist (override the general handler above)
 * -------------------------------------------------------------------------- */
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

/* --------------------------------------------------------------------------
 * 10. Forward declarations for OE types (SWIG-parsed, not compiled)
 * -------------------------------------------------------------------------- */
namespace OEChem {
    class OEMolBase;
    class OEMol;
}
namespace OEBio {
    class OEDesignUnit;
}

/* --------------------------------------------------------------------------
 * 11. Ignore problematic members before %include
 *
 *     StorageBackend.h pulls in <shared_mutex>, <unordered_map>, <thread>
 *     which SWIG cannot parse. We redeclare the classes manually instead.
 *
 *     Also ignore private clone constructors for metrics.
 * -------------------------------------------------------------------------- */

/* Suppress default constructors for types that don't have them */
%nodefaultctor OECluster::OEClusterError;
%nodefaultctor OECluster::MetricError;
%nodefaultctor OECluster::StorageError;
%nodefaultctor OECluster::DistanceMatrix;

/* Ignore Data() -- we expose _data_ptr() via %extend instead */
%ignore OECluster::StorageBackend::Data;
%ignore OECluster::DenseStorage::Data;
%ignore OECluster::MMapStorage::Data;
%ignore OECluster::SparseStorage::Data;

/* Ignore copy/move constructors and operators that SWIG cannot handle */
%ignore OECluster::MMapStorage::MMapStorage(const MMapStorage&);
%ignore OECluster::MMapStorage::operator=;

/* Ignore DistanceMetric::Clone -- returns unique_ptr, expose via %extend */
%ignore OECluster::DistanceMetric::Clone;

/* Ignore private constructors for metric clones */
%ignore OECluster::FingerprintMetric::FingerprintMetric(std::shared_ptr<const Impl>);

/* Ignore static mask parsing methods (not needed from Python) */
%ignore OECluster::FingerprintMetric::ParseAtomTypeMask;
%ignore OECluster::FingerprintMetric::ParseBondTypeMask;
%ignore OECluster::ROCSMetric::ROCSMetric(std::shared_ptr<const SharedData>, const Options&);
%ignore OECluster::SuperposeMetric::SuperposeMetric(std::shared_ptr<const SharedData>, const Options&);
/* Ignore SparseStorage internals that use unordered_map/shared_mutex/thread */
%ignore OECluster::SparseStorage::Entries;

/* --------------------------------------------------------------------------
 * 12. Redeclare StorageBackend hierarchy for SWIG
 *     (instead of %include "oecluster/StorageBackend.h" which has complex STL)
 * -------------------------------------------------------------------------- */
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

}  /* namespace OECluster */

/* --------------------------------------------------------------------------
 * 13. %extend blocks for StorageBackend and SparseStorage
 * -------------------------------------------------------------------------- */
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

/* --------------------------------------------------------------------------
 * 14. Error hierarchy
 * -------------------------------------------------------------------------- */
%include "oecluster/Error.h"

/* --------------------------------------------------------------------------
 * 15. DistanceMetric base class
 * -------------------------------------------------------------------------- */
%include "oecluster/DistanceMetric.h"

/* --------------------------------------------------------------------------
 * 16. ThreadPool (PIMPL -- only expose public interface)
 *     Redeclare to avoid pulling in complex threading headers.
 * -------------------------------------------------------------------------- */
namespace OECluster {

class ThreadPool {
public:
    explicit ThreadPool(size_t num_threads = 0);
    ~ThreadPool();
    size_t NumThreads() const;
    void Cancel();
    bool IsCancelled() const;
};

}  /* namespace OECluster */

/* --------------------------------------------------------------------------
 * 17. PDistOptions and pdist function
 *     Redeclare PDistOptions to handle the std::function member.
 * -------------------------------------------------------------------------- */
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

}  /* namespace OECluster */

/* --------------------------------------------------------------------------
 * 18. DistanceMatrix
 *     Redeclare to avoid unique_ptr constructor issues.
 * -------------------------------------------------------------------------- */
namespace OECluster {

class DistanceMatrix {
public:
    StorageBackend& Storage();
    const std::string& MetricName() const;
    const std::vector<std::string>& Labels() const;
    size_t NumItems() const;
    size_t NumPairs() const;
};

}  /* namespace OECluster */

/* --------------------------------------------------------------------------
 * 19. Metric classes
 * -------------------------------------------------------------------------- */
%include "oecluster/metrics/FingerprintMetric.h"
%include "oecluster/metrics/ROCSMetric.h"
%include "oecluster/metrics/SuperposeMetric.h"
