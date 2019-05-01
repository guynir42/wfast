function AT_CheckError( code )
% Converts an error code returned by the SDK into an error message and halts script execution.  
% AT_CheckError(code)
% Arguments:
%      code: the return code from the SDK
% Returns:
%      None
switch code
    case 0
        return %success!
    case 1
        error('AT_ERR_NONINITIALISED (1) Function called with an uninitialised handle');
    case 2
        error('AT_ERR_NOTIMPLEMENTED (2) Function not implemented for the chosen camera');
    case 3
        error('AT_ERR_READONLY (3) Feature is read only');
    case 4
        error('AT_ERR_NOTREADABLE (4) Feature is currently not readable');
    case 5
        error('AT_ERR_NOTWRITABLE (5) Feature is currently not writable');
    case 6
        error('AT_ERR_OUTOFRANGE (6) Value is outside the maximum and minimum limits');
    case 7
        error('AT_ERR_INDEXNOTAVAILABLE (7) Index is currently not available');
    case 8
        error('AT_ERR_INDEXNOTIMPLEMENTED (8) Index is not implemented for the chosen camera');
    case 9
        error('AT_ERR_EXCEEDEDMAXSTRINGLENGTH (9) Exceeded Maximum String Length');
    case 10
        error('AT_ERR_CONNECTION (10) Error connecting to or disconnecting from hardware');
    case 11
        error('AT_ERR_NODATA (11) No Data');
    case 12
        error('AT_ERR_INVALIDHANDLE (12) Invalid Handle');
    case 13
        error('AT_ERR_TIMEDOUT (13) The AT_WaitBuffer function timed out while waiting for data arrive in output queue');
    case 14
        error('AT_ERR_BUFFERFULL (14) The input queue has reached its capacity');
    case 15 
        error('AT_ERR_INVALIDSIZE (15) Invalid Size');
    case 16 
        error('AT_ERR_INVALIDALIGNMENT (16) Invalid Alignment');
    case 17
        error('AT_ERR_COMM (17) An error has occurred while communicating with hardware');
    case 18
        error('AT_ERR_STRINGNOTAVAILABLE (18) Index / String is not available');
    case 19
        error('AT_ERR_STRINGNOTIMPLEMENTED (19) Index / String is not implemented for the chosen camera');
    case 20
        error('AT_ERR_NULL_FEATURE (20) Null Feature');
    case 21
        error('AT_ERR_NULL_HANDLE (21) Null device handle passed to function');
    case 22
        error('AT_ERR_NULL_IMPLEMENTED_VAR (22) Null Implemented Var');
    case 23
        error('AT_ERR_NULL_READABLE_VAR (23) Null Readable Var');
    case 24
        error('AT_ERR_NULL_READONLY_VAR (24) Variable is Read-only');
    case 25
        error('AT_ERR_NULL_WRITABLE_VAR (25) Null Writable Var');
    case 26
        error('AT_ERR_NULL_MINVALUE (26) Null Minimum Value');
    case 27
        error('AT_ERR_NULL_MAXVALUE (27) Null Maximum Value');
    case 28
        error('AT_ERR_NULL_VALUE (28) Null Value');
    case 29
        error('AT_ERR_NULL_STRING (29) Null String');
    case 30
        error('AT_ERR_NULL_COUNT_VAR (30) Null Count Var');
    case 31
        error('AT_ERR_NULL_ISAVAILABLE_VAR (31) Null Isavailable Var');
    case 32
        error('AT_ERR_NULL_MAXSTRINGLENGTH (32) Null Max String Length');
    case 33
        error('AT_ERR_NULL_EVCALLBACK (33) Null Event Callback');
    case 34
        error('AT_ERR_NULL_QUEUE_PTR (34) Null Pointer Queue');
    case 35
        error('AT_ERR_NULL_WAIT_PTR (35) Null Wait Pointer');
    case 36
        error('AT_ERR_NULL_PTRSIZE (36) Null Pointersize');1
    case 37
        error('AT_ERR_NOMEMORY (37) No Memory');
    case 38
        error('AT_ERR_DEVICEINUSE (38) Function failed to connect to a device because it is already being used');
    otherwise
        error('AT_UNKNOWN_ERROR (%d) Unknown Error',code); 
        
  AT_FinaliseLibrary();

end

