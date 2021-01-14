function compile

    if isempty(getenv('HDF5'))
        error('Please install HDF5 and set the environmental variable to the right place...');
    end

    filename = 'write.cpp';
    dirname = fullfile(getenv('WFAST'), '+file/+mex');
    
    % str = 'mex CXXFLAGS="$CXXFLAGS -std=c++11 -static"';
    str = 'mex CXXFLAGS="$CXXFLAGS -std=c++11 -fpermissive -w"';
    str = [str ' -I' getenv('HDF5') '/include'];
    
    if ispc
        str = [str ' ' getenv('HDF5') '/lib/libhdf5.lib'];
        str = [str ' ' getenv('HDF5') '/lib/szip.lib'];
        str = [str ' ' getenv('HDF5') '/lib/zlib.lib'];
    else
        str = [str ' ' getenv('HDF5') '/lib/libhdf5.so'];
        str = [str ' ' getenv('HDF5') '/lib/libsz.so'];
        str = [str ' ' getenv('HDF5') '/lib/libz.so'];
    end
    
    str = [str ' ' fullfile(dirname, filename) ' -outdir ' dirname];
    str = [str ' ' fullfile(dirname, 'MyMatrix.cpp')];
    str = [str ' ' fullfile(dirname, 'MyAttribute.cpp')];
    str = [str ' ' fullfile(dirname, 'SaveData.cpp')];
    str = [str ' ' fullfile(dirname, 'SaveDataHDF5.cpp')];
    % add other datat save types here...

    str
    eval(str);
