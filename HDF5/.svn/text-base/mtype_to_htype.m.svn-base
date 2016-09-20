function htype = mtype_to_htype(mtype)
% MTYPE_TO_HTYPE converts matlab data types into the hdf5 equivalent.
%       HTYPE = MTYPE_TO_HTYPE(MTYPE)


if strmatch('H5T',mtype)
    htype = mtype;
    return
end



switch mtype
    case 'single'
        htype = 'H5T_NATIVE_FLOAT';
    case 'double'
        htype = 'H5T_NATIVE_DOUBLE';
    case 'char'
        htype = 'H5T_NATIVE_CHAR';
    case 'logical'
        htype = 'H5T_NATIVE_HBOOL';
    case 'int16'
        htype = 'H5T_NATIVE_SHORT';
    case 'uint16'
        htype = 'H5T_NATIVE_USHORT';
    case 'int32'
        htype = 'H5T_NATIVE_INT';
    case 'uint32'
        htype = 'H5T_NATIVE_UINT';
    case 'int64'
        htype = 'H5T_NATIVE_LONG';
    case 'uint64'
        htype = 'H5T_NATIVE_ULONG';
    case 'uint8'
        htype = 'H5T_STD_U8LE';
    otherwise
        error('Unrecognized type');
end
  