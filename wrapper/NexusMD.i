%module NexusMD

%{
#include "TaskDispatcher.h"
#include "Coords3D.h"
#include "SystemXMLParser.h"
%}

// Enable C++ standard library support
%include "std_string.i"
%include "std_shared_ptr.i"

// Expose shared_ptr for TaskDispatcher
%shared_ptr(TaskDispatcher)

// Include the classes
%include "TaskDispatcher.h"
%include "Coords3D.h"
%include "SystemXMLParser.h"
