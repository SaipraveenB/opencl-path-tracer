#include <renderer/camera.h>
// model camera written for GLFW.
// Automatically captures input from a GLFW service.
    // Two normalized vectors define direction.
    AbstractCamera::AbstractCamera( GLFWwindow* window ): glWindow(window){
    }
    ModelCamera::ModelCamera( GLFWwindow* window ):AbstractCamera( window ){
        //glfwSetInputMode( window, GLFW_CURSOR ,GLFW_CURSOR_DISABLED );
        iState = CAM_INACTIVE;
    }

    void ModelCamera::setSpeedX( float speed ){
        fSpeedX = speed;
    }

    void ModelCamera::setSpeedY( float speed ){
        fSpeedY = speed;
    }

    float ModelCamera::getSpeedX(){
        return fSpeedX;
    }

    float ModelCamera::getSpeedY(){
        return fSpeedY;
    }

    void ModelCamera::setRadius( float f ){
        fRadius = f;
    }
    void ModelCamera::setScrollSpeed( float f ){
        fScrollSpeed = f;
    }


    void AbstractCamera::setOrientation( glm::vec3 vec ){
        vOrientation = vec;
    }

    glm::mat4 AbstractCamera::getViewMatrix(){
        return glm::lookAt( vPos, vLookAt, vUp );
    }

    void ModelCamera::reset( glm::vec3 vDir ){
        vPos = glm::normalize( vDir ) * fRadius;

        vLookAt = glm::normalize( -vPos );

        // Do norm of lookAt x ( lookAt x orientation_base_vector ) to get the new Up vector.
        vUp = glm::normalize( glm::cross( vLookAt, glm::cross( vLookAt, vOrientation ) ) );

        //mView = glm::lookAt( vPos, vLookAt, vUp );
    }

    glm::vec3 ModelCamera::getEye(){
      return vPos;
    }
    glm::vec3 ModelCamera::getLookAt(){
      return vLookAt;
    }
    
    int ModelCamera::isChanged(){
        return (iState == CAM_ACTIVE);
        //return changed;
    }

    void ModelCamera::glfwHandleCursor( float fTimeDelta ){


      if( scrollState == SCROLL_INACTIVE ){
        if( glfwGetMouseButton( glWindow, GLFW_MOUSE_BUTTON_MIDDLE ) == GLFW_PRESS ){
          glfwSetInputMode( glWindow, GLFW_CURSOR ,GLFW_CURSOR_DISABLED );
          scrollState = SCROLL_ACTIVE;
          glfwSetCursorPos( glWindow, 0, 0 );
        }

      }else if( scrollState == SCROLL_ACTIVE ){
        if( glfwGetMouseButton( glWindow, GLFW_MOUSE_BUTTON_MIDDLE ) == GLFW_RELEASE ){
          glfwSetInputMode( glWindow, GLFW_CURSOR ,GLFW_CURSOR_NORMAL );
          scrollState = SCROLL_INACTIVE;

        }else{
          double x, y;
          glfwGetCursorPos( glWindow, &x, &y );
          glfwSetCursorPos( glWindow, 0, 0 );
          fRadius -= y*fScrollSpeed;
          reset( getEye() );
        }
      }

        if( iState == CAM_INACTIVE ){
          if( glfwGetMouseButton( glWindow, GLFW_MOUSE_BUTTON_LEFT ) == GLFW_PRESS ){
            glfwSetInputMode( glWindow, GLFW_CURSOR ,GLFW_CURSOR_DISABLED );
            iState = CAM_ACTIVE;
            glfwSetCursorPos( glWindow, 0, 0 );
          }
          return;
        }else if( iState == CAM_ACTIVE ){
          if( glfwGetMouseButton( glWindow, GLFW_MOUSE_BUTTON_LEFT ) == GLFW_RELEASE ){
            glfwSetInputMode( glWindow, GLFW_CURSOR ,GLFW_CURSOR_NORMAL );
            iState = CAM_INACTIVE;
            return;
          }
        }


        double x,y;
        glfwGetCursorPos( glWindow, &x, &y );
        glfwSetCursorPos( glWindow, 0, 0 );

        // Scale components to obtain degree space projections.
        double degX = fSpeedX * x;
        double degY = fSpeedY * y;

        // vUp is used as axis for rotating the point about virtual X
        // vUp x vLookAt isusedas pivot axis for virtual Y.

        // First compute the new differentially modified point.
        glm::vec3 newPos = vPos;
        newPos = glm::rotate( newPos, static_cast<float>(degX), vUp );
        newPos = glm::rotate( newPos, static_cast<float>(degY), glm::cross( vUp, vLookAt ) );
        vPos = newPos;

        // Now compute vUp and vLookAt.

        // just the inverse of the current position vector.
        vLookAt = glm::normalize( -vPos );

        // Do norm of lookAt x ( lookAt x orientation_base_vector ) to get the new Up vector.
        vUp = glm::normalize( glm::cross( vLookAt, glm::cross( vLookAt, vOrientation ) ) );

        //mView = glm::lookAt( vPos, vLookAt, vUp );
    }

CLCamera* ModelCamera::getCLCamera(){
    CLCamera* cam = new CLCamera();
    cam->vLookAt = glm::vec4( vLookAt, 0.0f );
    cam->vUp = glm::vec4( vUp, 0.0f );
    cam->vPos = glm::vec4( vPos, 0.0f );
    return cam;
}
