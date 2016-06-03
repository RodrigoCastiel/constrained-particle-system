#include "cps_program.h"

CPSProgram::CPSProgram() : GlutProgram()
{
  mScene = new CPSScene();
  mVideoRecorder = new VideoRecorder();
}

void CPSProgram::Init(int* argc, char* argv[], const char* windowTitle)
{
  GlutProgram::Init(argc, argv, windowTitle);
  GlutProgram::LoadShaders("../support");
  CPSProgram::InitScene(*argc, argv);
}

void CPSProgram::InitScene(int argc, char *argv[])
{
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT,  GL_NICEST);
  
  mScene->Init(mPipelineProgram, mProgramHandle);

  // TODO: initialize particle system.

}

// GLUT Callback methods ----------------------------------------------------------------

void CPSProgram::DisplayFunc()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  mScene->Render();
  glutSwapBuffers();
  mVideoRecorder->Update();
}

void CPSProgram::IdleFunc()
{
  mScene->Animate();
  glutPostRedisplay();
}

void CPSProgram::ReshapeFunc(int w, int h)
{
  GlutProgram::ReshapeFunc(w, h);

  mScene->ReshapeScreen(w, h);
  mVideoRecorder->UpdateSize(w, h);
}

void CPSProgram::PassiveMotionFunc(int x, int y)
{
  GlutProgram::PassiveMotionFunc(x, y);
}

void CPSProgram::MouseFunc(int button, int state, int x, int y)
{
  GlutProgram::MouseFunc(button, state, x, y);

  // NOTE: The following code was provided by the starter code.
  // keep track of whether CTRL and SHIFT keys are pressed
  switch (glutGetModifiers())
  {
    case GLUT_ACTIVE_ALT:
        mControlState = kEDIT;

        if (mMouse.mLftButton) // Clicking with left button in EDIT mode
        {
          mScene->OnMouseLeftClick(x, y, mWindowWidth, mWindowHeight);
        }
        else if (mMouse.mRgtButton) // Clicking with left button in EDIT mode
        {
          mScene->OnMouseRightClick(x, y, mWindowWidth, mWindowHeight);
        }

    break;

    case GLUT_ACTIVE_CTRL:
      mControlState = kTRANSLATE;
    break;

    case GLUT_ACTIVE_SHIFT:
      mControlState = kSCALE;
    break;

    // if CTRL and SHIFT are not pressed, we are in rotate mode
    default:
      mControlState = kROTATE;
    break;
  }
  // END.

}

// Drag function.
void CPSProgram::MotionFunc(int x, int y)
{
  Camera* camera = mScene->GetCurrentCamera();

  int mousePosDelta[2] = { x - mMouse.mPos[0], y - mMouse.mPos[1] };

  switch (mControlState)
  {
    // deform landscape
    case kEDIT:
      break;

    // translate the landscape
    case kTRANSLATE:
      if (mMouse.mLftButton)
      {
        camera->Translate(-mousePosDelta[0]/40.0f, +mousePosDelta[1]/10.0f, 0.0f);
      }
      if (mMouse.mRgtButton)
      {
        // control z translation via the right mouse button
        camera->Translate(0.0f, 0.0f, +mousePosDelta[1]/10.0f);
      }
      break;

    // rotate the landscape
    case kROTATE:
      if (mMouse.mLftButton)
      {
        camera->Rotate(-mousePosDelta[1]/100.0f, -mousePosDelta[0]/100.0f, 0.0f);
      }
      if (mMouse.mRgtButton)
      {
        // control z rotation via the right mouse button
        camera->Rotate(0.0f, 0.0f, mousePosDelta[1]/100.0f);
      }
      break;

    // scale the landscape
    case kSCALE:
      if (mMouse.mLftButton)
      {
        camera->Scale(+mousePosDelta[0]/100.0f, -mousePosDelta[1]/100.0f, 0.0f);
      }
      if (mMouse.mRgtButton)
      {
        // control z scaling via the right mouse button
        camera->Scale(0.0f, 0.0f, -mousePosDelta[1]/100.0f);
      }
      break;
  }
  // END.

  GlutProgram::MotionFunc(x, y);
}

void CPSProgram::KeyboardFunc(unsigned char key, int x, int y)
{
  GlutProgram::KeyboardFunc(key, x, y);

  switch (key)
  {
    case ' ':
      std::cout << "You pressed the spacebar." << std::endl;
      mVideoRecorder->ToggleRecord();
    break;

    case 'c':
      mScene->ChangeCamera();
    break;

    case 'f':
      glutFullScreen();
    break;

    case 'x':
      // take a screenshot
      mVideoRecorder->TakeScreenshot();
    break;
  }
}

// GLUT Callback methods end ------------------------------------------------------------