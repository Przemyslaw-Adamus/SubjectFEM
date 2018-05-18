#include "Element.h"

Element::Element()
{
	this->elementId = 0;
	for (int i = 0;i < 4;i++) 
	{
		this->nodeId[i] = 0;
	}
	for (int i = 0;i < 4;i++)
	{
		boundaryConditions[i] = false;
	}
}

