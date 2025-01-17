%module Nexus

%{
#include "TaskDispatcher.h"
%}

%include "std_string.i"      // Enable std::string support
%include "std_shared_ptr.i"  // Enable shared_ptr support
%shared_ptr(TaskDispatcher)  // Use shared_ptr for TaskDispatcher

%include "TaskDispatcher.h"

