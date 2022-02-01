# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c

# Include any dependencies generated for this target.
include src/CMakeFiles/NucleiSegregation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/NucleiSegregation.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/NucleiSegregation.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/NucleiSegregation.dir/flags.make

src/CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o: src/CMakeFiles/NucleiSegregation.dir/flags.make
src/CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o: src/nuclei_segregation.c
src/CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o: src/CMakeFiles/NucleiSegregation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT src/CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o -MF CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o.d -o CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o -c /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/nuclei_segregation.c

src/CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.i"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/nuclei_segregation.c > CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.i

src/CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.s"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/nuclei_segregation.c -o CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.s

src/CMakeFiles/NucleiSegregation.dir/params.c.o: src/CMakeFiles/NucleiSegregation.dir/flags.make
src/CMakeFiles/NucleiSegregation.dir/params.c.o: src/params.c
src/CMakeFiles/NucleiSegregation.dir/params.c.o: src/CMakeFiles/NucleiSegregation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/CMakeFiles/NucleiSegregation.dir/params.c.o"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT src/CMakeFiles/NucleiSegregation.dir/params.c.o -MF CMakeFiles/NucleiSegregation.dir/params.c.o.d -o CMakeFiles/NucleiSegregation.dir/params.c.o -c /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/params.c

src/CMakeFiles/NucleiSegregation.dir/params.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/NucleiSegregation.dir/params.c.i"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/params.c > CMakeFiles/NucleiSegregation.dir/params.c.i

src/CMakeFiles/NucleiSegregation.dir/params.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/NucleiSegregation.dir/params.c.s"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/params.c -o CMakeFiles/NucleiSegregation.dir/params.c.s

src/CMakeFiles/NucleiSegregation.dir/output.c.o: src/CMakeFiles/NucleiSegregation.dir/flags.make
src/CMakeFiles/NucleiSegregation.dir/output.c.o: src/output.c
src/CMakeFiles/NucleiSegregation.dir/output.c.o: src/CMakeFiles/NucleiSegregation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/CMakeFiles/NucleiSegregation.dir/output.c.o"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT src/CMakeFiles/NucleiSegregation.dir/output.c.o -MF CMakeFiles/NucleiSegregation.dir/output.c.o.d -o CMakeFiles/NucleiSegregation.dir/output.c.o -c /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/output.c

src/CMakeFiles/NucleiSegregation.dir/output.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/NucleiSegregation.dir/output.c.i"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/output.c > CMakeFiles/NucleiSegregation.dir/output.c.i

src/CMakeFiles/NucleiSegregation.dir/output.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/NucleiSegregation.dir/output.c.s"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/output.c -o CMakeFiles/NucleiSegregation.dir/output.c.s

src/CMakeFiles/NucleiSegregation.dir/display_image.cpp.o: src/CMakeFiles/NucleiSegregation.dir/flags.make
src/CMakeFiles/NucleiSegregation.dir/display_image.cpp.o: src/display_image.cpp
src/CMakeFiles/NucleiSegregation.dir/display_image.cpp.o: src/CMakeFiles/NucleiSegregation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/NucleiSegregation.dir/display_image.cpp.o"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/NucleiSegregation.dir/display_image.cpp.o -MF CMakeFiles/NucleiSegregation.dir/display_image.cpp.o.d -o CMakeFiles/NucleiSegregation.dir/display_image.cpp.o -c /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/display_image.cpp

src/CMakeFiles/NucleiSegregation.dir/display_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NucleiSegregation.dir/display_image.cpp.i"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/display_image.cpp > CMakeFiles/NucleiSegregation.dir/display_image.cpp.i

src/CMakeFiles/NucleiSegregation.dir/display_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NucleiSegregation.dir/display_image.cpp.s"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/display_image.cpp -o CMakeFiles/NucleiSegregation.dir/display_image.cpp.s

# Object files for target NucleiSegregation
NucleiSegregation_OBJECTS = \
"CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o" \
"CMakeFiles/NucleiSegregation.dir/params.c.o" \
"CMakeFiles/NucleiSegregation.dir/output.c.o" \
"CMakeFiles/NucleiSegregation.dir/display_image.cpp.o"

# External object files for target NucleiSegregation
NucleiSegregation_EXTERNAL_OBJECTS =

src/NucleiSegregation: src/CMakeFiles/NucleiSegregation.dir/nuclei_segregation.c.o
src/NucleiSegregation: src/CMakeFiles/NucleiSegregation.dir/params.c.o
src/NucleiSegregation: src/CMakeFiles/NucleiSegregation.dir/output.c.o
src/NucleiSegregation: src/CMakeFiles/NucleiSegregation.dir/display_image.cpp.o
src/NucleiSegregation: src/CMakeFiles/NucleiSegregation.dir/build.make
src/NucleiSegregation: /usr/lib/libopencv_gapi.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_stitching.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_alphamat.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_aruco.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_barcode.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_bgsegm.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_bioinspired.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_ccalib.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_cvv.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_dnn_objdetect.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_dnn_superres.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_dpm.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_face.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_freetype.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_fuzzy.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_hdf.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_hfs.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_img_hash.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_intensity_transform.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_line_descriptor.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_mcc.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_quality.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_rapid.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_reg.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_rgbd.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_saliency.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_stereo.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_structured_light.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_superres.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_surface_matching.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_tracking.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_videostab.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_viz.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_wechat_qrcode.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_xfeatures2d.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_xobjdetect.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_xphoto.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_shape.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_highgui.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_datasets.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_plot.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_text.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_ml.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_phase_unwrapping.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_optflow.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_ximgproc.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_video.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_videoio.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_imgcodecs.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_objdetect.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_calib3d.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_dnn.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_features2d.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_flann.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_photo.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_imgproc.so.4.5.5
src/NucleiSegregation: /usr/lib/libopencv_core.so.4.5.5
src/NucleiSegregation: src/CMakeFiles/NucleiSegregation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable NucleiSegregation"
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NucleiSegregation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/NucleiSegregation.dir/build: src/NucleiSegregation
.PHONY : src/CMakeFiles/NucleiSegregation.dir/build

src/CMakeFiles/NucleiSegregation.dir/clean:
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src && $(CMAKE_COMMAND) -P CMakeFiles/NucleiSegregation.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/NucleiSegregation.dir/clean

src/CMakeFiles/NucleiSegregation.dir/depend:
	cd /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src /data/Perkins/TongueSTOmics/scripts/subsetting/nuclei_segregation_c/src/CMakeFiles/NucleiSegregation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/NucleiSegregation.dir/depend

