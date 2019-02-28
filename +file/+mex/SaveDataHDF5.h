#ifndef SAVEDATAHDF5_H
#define SAVEDATAHDF5_H

#include "SaveData.h"
#include "hdf5.h"

#define DESC_LENGTH 256

class SaveDataHDF5: public SaveData {
	
public:

	class MyFilePointer{// keeps a file id and closes it when destroyed
	
	public: 
	
		hid_t id;
	
		MyFilePointer(const char *filename){ id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); }
		~MyFilePointer(){ H5Fclose(id); }
	
	};

	class MyDataspace{// keeps a dataspace id and closes it when destroyed
	
	public: 
		
		hid_t id;
		hid_t data_type;
		
		MyDataspace(MyMatrix matrix);
		MyDataspace(MyAttribute attribute);
		~MyDataspace(){ H5Tclose(data_type); H5Sclose(id); }
		
	};

	class MyDataset{// keeps a dataset id and closes it when destroyed, also takes care of sizes and deflate
		
	public: 
		
		hid_t id;
		
		MyDataset(MyFilePointer &file, const char *location, MyMatrix matrix, MyDataspace &dataspace);
		~MyDataset(){ H5Dclose(id); }
		
	};
	
	class MyGroup{// keeps a group_id and closes it when destroyed
	
	public:
	
		hid_t id;
		
		MyGroup(MyFilePointer &file, const char *location);
		~MyGroup(){ H5Gclose(id); }
	
	};
	
	MyFilePointer file;
	
	SaveDataHDF5(const char *filename);
	~SaveDataHDF5();
	virtual void writeData();
	virtual void writeMatrix(MyFilePointer &file, MyMatrix matrix);
	virtual void writePars();
	virtual void writeAttribute(MyAttribute att, hid_t group_or_dataset_id);


	
};

#endif