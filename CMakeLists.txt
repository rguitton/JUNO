# CMAKE_BINARY_DIR est le répertoire d’où sera lancée la commande cmake, à savoir build
# S’assure que la version de make est suffisamment récente
cmake_minimum_required(VERSION 2.6)
# Indique le nom du projet
project(demo_cpp)
# Impose le répertoire de création du fichier exécutable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
# Impose le répertoire des fichiers source
set(CMAKE_CURRENT_SOURCE_DIR ${CMAKE_BINARY_DIR}/../src)
# Ajoute le répertoire include des fichiers .h au chemin de recherche des #include "xx.h" 
include_directories(${CMAKE_BINARY_DIR}/../include)

# ************************ Start: Specific to ROOT codes *************************

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}/cmake")

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT COMPONENTS RIO Net)


# ************************ End: Specific to ROOT codes *************************


# Impose certaines options de compilation
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra -Wpedantic $(root-config --cflags --libs)")
# Génère la liste des fichiers .cpp devant être compilés
file(
   GLOB_RECURSE
   srcFiles
   ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)
# Crée l'exécutable monExe à partir de la liste des fichiers .cpp
add_executable(
   monExe
   ${srcFiles}
)

