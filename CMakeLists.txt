# Copyright (c) 2021 Idiap Research Institute, http://www.idiap.ch/
# Written by Rudolf A. Braun <rbraun@idiap.ch>
# 
# This file is part of icassp-oov-recognition
# 
# icassp-oov-recognition is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as
# published by the Free Software Foundation.
# 
# icassp-oov-recognition is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with icassp-oov-recognition. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.9.5)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "-Ofast ")

include_directories("/path/to/openfst-1.6.7/include")
include_directories("libs/")

add_subdirectory(libs/pybind11)
pybind11_add_module(fast libs/fast.cc libs/fst-wrapper.cc)

set_target_properties(fast PROPERTIES LIBRARY_OUTPUT_NAME "fast")

target_link_libraries(fast PRIVATE "-L/path/to/openfst-1.6.7/lib" -lfstscript -lfst)
