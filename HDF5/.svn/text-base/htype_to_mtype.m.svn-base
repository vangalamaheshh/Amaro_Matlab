function mtype = htype_to_mtype(htype)
% HTYPE_TO_DTYPE converts hdf5 data type into the matlab equivalent.
%       MTYPE = HTYPE_TO_MTYPE(HTYPE)
%
%   For unexpected errors using this function, add the hdf5 to matlab
%   conversion to the switch statement.
%
%       Revisions:
%           22 Nov 07  -- Function added by Jen Dobson
%

switch htype
    case {'H5T_NATIVE_FLOAT','H5T_IEEE_F32LE'}
        mtype = 'single';
    case {'H5T_NATIVE_DOUBLE','H5T_IEEE_F64LE'}
        mtype = 'double';
    case 'H5T_NATIVE_CHAR'
        mtype = 'char';
    case 'H5T_NATIVE_HBOOL'
        mtype = 'logical';
    case'H5T_NATIVE_SHORT';
        mtype = 'int16';
    case  'H5T_NATIVE_USHORT';
        mtype = 'uint16';
    case  'H5T_NATIVE_INT';
        mtype = 'int32';
    case  'H5T_NATIVE_UINT';
        mtype = 'uint32';
    case  'H5T_NATIVE_LONG';
        mtype = 'int64';
    case  'H5T_NATIVE_ULONG';
        mtype = 'uint64';
    case 'H5T_STD_U8LE';
        mtype = 'uint8';

    otherwise
        error('Unrecognized Data Type (consider adding it to htype_to_mtype.m!)')
end
