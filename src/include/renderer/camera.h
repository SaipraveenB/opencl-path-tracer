#include "stdgl.h"


class AbstractCamera{
protected:
  glm::vec3 vLookAt;
  // an additional two components... note that not only is this normalized but it's also constrained to be perpendicular to lookAt. hence only one degree of freedom is offered by this component
  glm::vec3 vUp;

  // describes spatial point position of the camera.
  glm::vec3 vPos;

  glm::mat4 mView;

  glm::vec3 vOrientation;

  GLFWwindow* glWindow;

  int iState;
  int scrollState;

  static const int CAM_ACTIVE = 1;
  static const int CAM_INACTIVE = 0;

  static const int SCROLL_ACTIVE = 1;
  static const int SCROLL_INACTIVE = 0;

public:
  AbstractCamera( GLFWwindow* );
  glm::mat4 getViewMatrix();
  virtual void glfwHandleCursor( float ) = 0;
  void setOrientation( glm::vec3 );


};

struct CLCamera{
    glm::vec4 vPos;
    glm::vec4 vLookAt;
    glm::vec4 vUp;
};

class ModelCamera : public AbstractCamera{
    // Two normalized vectors define direction.
    // Three component vector..  but it's normalized hence only two degrees of freedom
    /*glm::vec3 vLookAt;
    // an additional two components... note that not only is this normalized but it's also constrained to be perpendicular to lookAt. hence only one degree of freedom is offered by this component
    glm::vec3 vUp;

    // describes spatial point position of the camera.
    glm::vec3 vPos;

    glm::mat4 mView;*/
    // Additional control variables..
    // Handle rotation speed in degrees of rotation/screen space pixel.
    float fSpeedX;
    float fSpeedY;
    float fScrollSpeed;
    // Clamp to a radius around the center of the world.
    float fRadius;
    // Single vector denoting the line towards which the camera is to be oriented.
    // Used to recalculate vUp vector;
    //glm::vec3 vOrientation;

    // UI Binder variables.
    //GLFWwindow* glWindow;

    public:
    ModelCamera( GLFWwindow* );
    void setSpeedX( float );
    void setSpeedY( float );
    void setScrollSpeed( float );
    float getSpeedX();
    float getSpeedY();
    void setRadius( float );
    //void setOrientation( glm::vec3 );
    void reset( glm::vec3 );
    //glm::mat4 getViewMatrix();
    void glfwHandleCursor( float );
    glm::vec3 getEye();
    glm::vec3 getLookAt();
    int isChanged();
    CLCamera* getCLCamera();
};

class FirstPersonCamera : public AbstractCamera{
  float fOmegaX;
  float fOmegaY;
  // Clamp to a radius around the center of the world.
  float fSpeed;
  // Single vector denoting the line towards which the camera is to be oriented.
  // Used to recalculate vUp vector;
  //glm::vec3 vOrientation;

  // UI Binder variables.
  //GLFWwindow* glWindow;

public:
  FirstPersonCamera( GLFWwindow* );
  void setOmegaX( float );
  void setOmegaY( float );
  void setSpeed( float );
  //void setOrientation( glm::vec3 );
  void reset( glm::vec3 );
  //glm::mat4 getViewMatrix();
  void glfwHandleCursor( float );
};
