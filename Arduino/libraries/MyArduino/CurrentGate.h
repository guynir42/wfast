#include <Arduino.h>

class CurrentGate {

public:

	CurrentGate(short int pin1, short int pin2);
	
	bool getState() const;
	bool getStatus() const;
	bool getInverse() const;
	
	void setState(bool state);
	void setInverse(bool inverse=1);
	

protected:

	void allowCurrent();
	void blockCurrent();

	short int _pin1;
	short int _pin2;
	
	bool _state=0;
	bool _inverse=0;
	
};
