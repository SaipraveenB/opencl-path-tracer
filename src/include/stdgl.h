#ifndef GL_MAIN
#define GL_MAIN

//#include <GL/glew.h> // include GLEW and new version of GL on Windows
#ifdef __linux__
#include <GL/glew.h>
#include <GL/glext.h>
#endif

#ifdef __APPLE__
#include <GL/glew.h>
#endif

#include <GLFW/glfw3.h> // GLFW helper library
#include <glm/glm.hpp>

// glm::translate, glm::rotate, glm::scale
#include <glm/gtc/matrix_transform.hpp>
// glm::value_ptr
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#endif
