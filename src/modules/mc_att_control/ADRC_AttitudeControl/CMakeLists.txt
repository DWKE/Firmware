px4_add_library(ADRC_AttitudeControl
        ADRC_AttitudeControl.cpp
)
target_include_directories(ADRC_AttitudeControl
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
)

px4_add_unit_gtest(SRC ADRC_AttitudeControlTest.cpp LINKLIBS ADRC_AttitudeControl)
