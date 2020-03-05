#include "SaveDataHDF5.h"

SaveDataHDF5::SaveDataHDF5(const char *filename) : SaveData(filename), file(filename){
		
}

SaveDataHDF5::~SaveDataHDF5(){
	
	
}

void SaveDataHDF5::writeData(){ // write to file all the main data sets (images, timestamps, etc...)
	
	// MyFilePointer file(filename);
	
	// mexPrintf("writing data into HDF5\n");
	
	// if any of these is empty, they will get skipped...
	writeMatrix(file, images);
	writeMatrix(file, timestamps);
	writeMatrix(file, cutouts);
	writeMatrix(file, positions);
	writeMatrix(file, coordinates);
	writeMatrix(file, magnitudes);
	writeMatrix(file, temperatures);
	writeMatrix(file, cutouts_bg);
	writeMatrix(file, positions_bg);
	writeMatrix(file, stack);
	writeMatrix(file, psfs);
	
	writeMatrix(file, fluxes);
	writeMatrix(file, errors);
	writeMatrix(file, areas);
	writeMatrix(file, backgrounds);
	writeMatrix(file, variances);
	writeMatrix(file, offsets_x);
	writeMatrix(file, offsets_y);
	writeMatrix(file, centroids_x);
	writeMatrix(file, centroids_y);
	writeMatrix(file, widths);
	writeMatrix(file, bad_pixels);
	writeMatrix(file, flags);
	
}

void SaveDataHDF5::writeHeader(){ // write to file the head.Header object
	
	for(int i=0;i<parameter_addresses_vector.size();i++){

		MyGroup g(file, parameter_addresses_vector[i].c_str());
		
		for(int j=0;j<parameter_attributes_2D_vector[i].size();j++){
			
			writeAttribute(parameter_attributes_2D_vector[i][j], g.id);
			
		}
		
	}// for i
	
}

void SaveDataHDF5::writeMatrix(MyFilePointer &file, MyMatrix matrix){ // write a single dataset
	
	if(matrix.is_empty()) return;
	
	char location[STRLN];
	snprintf(location, STRLN, "/%s", matrix.data_name);
	
	MyDataspace dataspace(matrix);
	MyDataset dataset(file, location, matrix, dataspace);
	
	for(int i=0;i<matrix.attributes.size();i++) writeAttribute(matrix.attributes[i], dataset.id);
	
}

void SaveDataHDF5::writeAttribute(MyAttribute att, hid_t group_or_dataset_id){ // write a single attribute
	
	MyDataspace dataspace(att);
	
	hid_t id=H5Acreate2(group_or_dataset_id, att.att_name, dataspace.data_type, dataspace.id, H5P_DEFAULT, H5P_DEFAULT);

	int status=0;
	if(att.is_scalar) status=H5Awrite(id, dataspace.data_type, &att.scalar);
	if(att.is_vec) status=H5Awrite(id, dataspace.data_type, att.vec.data());
	if(att.is_str) status=H5Awrite(id, dataspace.data_type, att.str.c_str());
	
	if(status<0){ 
		mex_flag[2]=-1; 
		char str[10];
		if(att.is_empty()) strncpy(str, "empty", 10);
		else if(att.is_scalar) strncpy(str, "scalar", 10);
		else if(att.is_vec) strncpy(str, "vector", 10); 
		else if(att.is_str) strncpy(str, "string", 10);
		mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveDataHDF5:writeAttribute", "Something wrong when writing %s attribute '%s'", str, att.att_name); 
	}
	
	H5Tclose(dataspace.data_type);
	H5Aclose(id);
	
}

// definitions for MyDataset and MyDataspace and MyGroup
SaveDataHDF5::MyDataspace::MyDataspace(MyMatrix matrix){ 

	id=H5Screate_simple(matrix.ndims, matrix.dims_c, NULL); 
	data_type=H5T_NATIVE_DOUBLE;
	if(matrix.is_uint16()) data_type=H5T_NATIVE_UINT16;
	if(matrix.is_float()) data_type=H5T_NATIVE_FLOAT;
	
	if(id<0){ 
		mex_flag[2]=-1; 
		mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveDataHDF5:create_dataspace", "Something wrong when creating dataspace for '%s'", matrix.data_name); 
	}
	
}

SaveDataHDF5::MyDataspace::MyDataspace(MyAttribute attribute){
	
	// default is to use a scalar...
	// size_t dims_temp=0;
	data_type=H5T_NATIVE_DOUBLE;
	
	if(attribute.is_empty()) id=H5Screate(H5S_NULL);
	else if(attribute.is_scalar) id=H5Screate(H5S_SCALAR);
	else if(attribute.is_vec){ // allow for vector attributes
		id=H5Screate(H5S_SIMPLE); 
		size_t dims[1]={attribute.vec.size()};
		H5Sset_extent_simple(id, 1, dims, NULL);
	}
	else if(attribute.is_str){
		
		id=H5Screate(H5S_SIMPLE); 
		size_t dims[1]={1};
		H5Sset_extent_simple(id, 1, dims, NULL);
		data_type = H5Tcopy(H5T_C_S1);
		H5Tset_size(data_type, attribute.str.length());
		// H5Tset_size(data_type, H5T_VARIABLE);
		H5Tset_strpad(data_type, H5T_STR_NULLTERM);
		
	}
	
	// id=H5Screate_simple(1, &dims_temp, NULL);
	
	if(id<0){ 
		mex_flag[2]=-1; 
		mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveDataHDF5:create_dataspace", "Something wrong when creating dataspace for '%s'", attribute.att_name); 
	}
	
}

SaveDataHDF5::MyDataset::MyDataset(MyFilePointer &file, const char *location, MyMatrix matrix, MyDataspace &dataspace){
	
	size_t chunk=chunk_size; // verify the chunk size is not bigger than the matrix size
	if(chunk>matrix.cols) chunk=matrix.cols;
	if(chunk>matrix.rows) chunk=matrix.rows;
	
	int data_type=H5T_NATIVE_DOUBLE; // default data type is double
	if(matrix.is_uint16()) data_type=H5T_NATIVE_USHORT; 
	if(matrix.is_float()) data_type=H5T_NATIVE_FLOAT;
	
	hid_t plist_id=H5Pcreate(H5P_DATASET_CREATE);
	
	if (matrix.use_deflate && deflate>0){ // if this matrix needs to be deflated and if we are using deflate (in general)
		
		hsize_t chunk_dims_matlab[4]={chunk, chunk, 1, 1};
		
		if(chunk<chunk_size) chunk_dims_matlab[2] = matrix.frames; // for small datasets like cutouts we want all frames of the same cutout chunked together...
		
		hsize_t chunk_dims_c[4]={1};
		for(int i=0;i<matrix.ndims;i++) chunk_dims_c[i]=chunk_dims_matlab[matrix.ndims-1-i];
		
		int status=0;
		
		status = H5Pset_chunk (plist_id, matrix.ndims, chunk_dims_c);
		if(status<0){ 
			mex_flag[2]=-1; 
			mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:set_chunk", "Something wrong when setting chunk for %s!", matrix.data_name);
		}
		
		status = H5Pset_deflate (plist_id, deflate);
		if(status<0){ 
			mex_flag[2]=-1; 
			mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:set_deflate", "Something wrong when setting deflate! for %s", matrix.data_name); 
		}
		
	}
	
	id=H5Dcreate(file.id, location, data_type, dataspace.id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	
	// go on and write the matrix to file...
	int status=0;
	if(matrix.is_double()){// write double data
	
		status = H5Dwrite(id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix.matrix_double);
		if(debug_bit>2) mexPrintf("saving %s: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", matrix.data_name, file.id, dataspace.id, id, status);
		if(status<0){ 
			mex_flag[2]=-1; 
			mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveDataHDF5:write_dataset", "Something wrong when writing dataset for %s (double)!", matrix.data_name); 
		}
			
	}
	else if(matrix.is_uint16()){// write uint16 data
	
		status = H5Dwrite(id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix.matrix_uint16);
		if(debug_bit>2) mexPrintf("saving %s: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", matrix.data_name, file.id, dataspace.id, id, status);
		if(status<0){ 
			mex_flag[2]=-1; 
			mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveDataHDF5:write_dataset", "Something wrong when writing dataset for %s (uint16)!", matrix.data_name); 
		}
	}
	else if(matrix.is_float()){// write float data
	
		status = H5Dwrite(id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix.matrix_float);
		if(debug_bit>2) mexPrintf("saving %s: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", matrix.data_name, file.id, dataspace.id, id, status);
		if(status<0){ 
			mex_flag[2]=-1; 
			mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveDataHDF5:write_dataset", "Something wrong when writing dataset for %s (float)!", matrix.data_name); 
		}
	}
	
	
}

SaveDataHDF5::MyGroup::MyGroup(MyFilePointer &file, const char *location){
	
	id=H5Gcreate2(file.id, location, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
}




