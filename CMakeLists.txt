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

### ROOT configuration de "Joao Pedro Athayde Marcondes de Andre" <jpandre@iphc.cnrs.fr>; 
#find_package(ROOT QUIET COMPONENTS)
if(ROOT_FOUND)
     message("${Cyan}-- CMake configured ROOT found${ColourReset}")
     include(${ROOT_USE_FILE})
     include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
else(ROOT_FOUND)
     find_program ( ROOT_CONFIG_EXEC root-config )
     if( ROOT_CONFIG_EXEC )
         message("${Cyan}-- root-config found:
${ROOT_CONFIG_EXEC}${ColourReset}")
         set( ROOT_FOUND ON )
         set( ROOT_FOUND_ROOTCONFIG ON )
         execute_process(
             COMMAND sh "${ROOT_CONFIG_EXEC}" --incdir
             OUTPUT_VARIABLE ROOT_INCLUDE_DIRS
OUTPUT_STRIP_TRAILING_WHITESPACE
             )
         include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
         execute_process(
             COMMAND sh "${ROOT_CONFIG_EXEC}" --auxcflags
             OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE
             )
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

         execute_process(
             COMMAND sh "${ROOT_CONFIG_EXEC}" --libdir
             OUTPUT_VARIABLE ROOT_LIB_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
             )
         execute_process(
             COMMAND sh "${ROOT_CONFIG_EXEC}" --libs
             OUTPUT_VARIABLE ROOT_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE
             )
         STRING(REGEX REPLACE "-l([^ ]*)" "lib\\1.so" ROOT_LIBS
"${ROOT_LIBS}")
         STRING(REGEX MATCHALL "lib[^ ]*.so" ROOT_LIBS "${ROOT_LIBS}")
         SET(ROOT_LIBRARIES "")
         foreach(lib ${ROOT_LIBS})
             find_library(ROOT_${lib}_LIBRARY
                 NAMES ${lib}
                 PATHS ${ROOT_LIB_DIR}
                 )
             # message("${lib} => ${ROOT_${lib}_LIBRARY}")
             if(ROOT_${lib}_LIBRARY)
                 set(ROOT_LIBRARIES ${ROOT_LIBRARIES}
${ROOT_${lib}_LIBRARY})
             endif(ROOT_${lib}_LIBRARY)
         ENDforeach()
     else( ROOT_CONFIG_EXEC )
         message("${Red}-- root-config not found in PATH. Using
ROOTSYS-based configuration${ColourReset}")
         list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
         find_package(ROOT REQUIRED COMPONENTS)
         include(${ROOT_USE_FILE})
     endif( ROOT_CONFIG_EXEC )
endif(ROOT_FOUND)

# ************************ End: Specific to ROOT codes *************************

# Définit certaines options de compilation
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic ${ROOT_CXX_FLAGS}")

# Génère la liste des fichiers .cpp à compiler
file(GLOB_RECURSE srcFiles ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Crée l'exécutable à partir des fichiers .cpp
add_executable(monExe ${srcFiles})

# Liens avec les bibliothèques ROOT
target_link_libraries(monExe ${ROOT_LIBRARIES})
