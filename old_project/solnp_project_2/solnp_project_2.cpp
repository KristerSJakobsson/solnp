// solnp_project_2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "solnp.hpp"
#include "tests.cpp"


int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

