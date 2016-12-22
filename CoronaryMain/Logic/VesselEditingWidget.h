//#ifndef __VTK_WRAP__


#ifndef __VesselEditingWidget_h
#define __VesselEditingWidget_h


#include "QVTKWidget.h"


class QVesselEditingWidget : public QVTKWidget
{
	Q_OBJECT
public:
	QVesselEditingWidget(QVTKWidget *parent = 0, const char *name = 0){};
	~QVesselEditingWidget() {}

protected:
	virtual void showEvent(QShowEvent* e);

public:
};



#endif


//#endif //__VTK_WRAP__