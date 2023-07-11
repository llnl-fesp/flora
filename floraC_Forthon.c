/* Created by Ron Cohen for UEDGE */
/* This is a stub module which calls the init functions of all of            */
/* the modules that are part of Flora. This is needed since the modules      */
/* depend on each other and so must be incorporated into one shared          */
/* object file.                                                              */
#include <Python.h>
#define NPY_NO_DEPRECATED_API 8
#include <numpy/arrayobject.h>


#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

static PyObject *ErrorObject;

#if PY_MAJOR_VERSION < 3
  extern PyMODINIT_FUNC initglrpy(void);
  extern PyMODINIT_FUNC initutlpy(void);
#endif

/* ######################################################################### */
/* # Method list                                                             */
static struct PyMethodDef floraC_methods[] = {
  {NULL,NULL}};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "floraC", /* m_name */
  "floraC", /* m_doc */
  -1,                  /* m_size */
  floraC_methods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
  };
#endif



/* ######################################################################### */
/* # The initialization function                                             */
#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_floraC(void)
#else
PyMODINIT_FUNC initfloraC(void)
#endif
{

  PyObject *m, *d;
  /*PyObject *pystdout;*/
#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
#else
  m = Py_InitModule("floraC", floraC_methods);
#endif

  d = PyModule_GetDict(m);
  ErrorObject = PyErr_NewException("floraC.error",NULL,NULL);
  PyDict_SetItemString(d, "error", ErrorObject);
  if (PyErr_Occurred())
    Py_FatalError("can not initialize module floraC");

  /*pystdout = PySys_GetObject("stdout");*/
  /*PyFile_WriteString("Forthon edition\n",pystdout);*/

  import_array();

#if PY_MAJOR_VERSION < 3
  initglrpy();
  initutlpy();
#endif
#if PY_MAJOR_VERSION >= 3
  return m;
#else
  return;
#endif

}


