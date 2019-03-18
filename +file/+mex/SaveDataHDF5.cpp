#include "SaveDataHDF5.h"

SaveDataHDF5::SaveDataHDF5(const char *filename) : SaveData(filename), file(filename){
		
}

SaveDataHDF5::~SaveDataHDF5(){
	
	
}

void SaveDataHDF5::writeData(){ // write to file all the main data sets (images, timestamps, etc...)
	
	// MyFilePointer file(filename);
	
	// if any of these is empty, they will get skipped...
	writeMatrix(file, images);
	writeMatrix(file, cutouts);
	writeMatrix(file, positions);
	writeMatrix(file, stack);
	writeMatrix(file, timestamps);
	writeMatrix(file, psfs);
	writeMatrix(file, lightcurves);
	
}

void SaveDataHDF5::writePars(){ // write to file the head.Parameters object
	
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
		mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveDataHDF5:writeAttribute", "Something wrong when writing attribute '%s'", att.att_name); 
	}
	
	H5Tclose(dataspace.data_type);
	H5Aclose(id);
	
}

// definitions for MyDataset and MyDataspace and MyGroup
SaveDataHDF5::MyDataspace::MyDataspace(MyMatrix matrix){ 

	id=H5Screate_simple(matrix.ndims, matrix.dims_c, NULL); 
	data_type=H5T_NATIVE_DOUBLE;
	if(matrix.is_uint16()) data_type=H5T_NATIVE_UINT16;
	
	if(id<0){ 
		mex_flag[2]=-1; 
		mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveDataHDF5:create_dataspace", "Something wrong when creating dataspace for '%s'", matrix.data_name); 
	}
	
}

SaveDataHDF5::MyDataspace::MyDataspace(MyAttribute attribute){
	
	// default is to use a scalar...
	size_t dims_temp=1;
	data_type=H5T_NATIVE_DOUBLE;
	
	if(attribute.is_vec) dims_temp=attribute.vec.size(); // allow for vector attributes
	
	id=H5Screate_simple(1, &dims_temp, NULL);
	
	if(attribute.is_str){
		
		data_type = H5Tcopy(H5T_C_S1);
		H5Tset_size(data_type, attribute.str.length());
		// H5Tset_size(data_type, H5T_VARIABLE);
		H5Tset_strpad(data_type, H5T_STR_NULLTERM);
		
	}
	
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
	
	hid_t plist_id=H5Pcreate(H5P_DATASET_CREATE);
	
	if (matrix.use_deflate && deflate>0){ // if this matrix needs to be deflated and if we are using deflate (in general)
		
		hsize_t chunk_dims_c[4]={1};
		for(int i=0;i<matrix.ndims;i++) chunk_dims_c[matrix.ndims-1-i]=chunk_size;
		
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
	
	
}

SaveDataHDF5::MyGroup::MyGroup(MyFilePointer &file, const char *location){
	
	id=H5Gcreate2(file.id, location, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
}




