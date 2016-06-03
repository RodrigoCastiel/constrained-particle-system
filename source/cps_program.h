/******************************************+
*                                          *
*  Author: Rodrigo Castiel                 *
*  Email: rcrs2@cin.ufpe.br                *
*         castielr@usc.edu                 *
*                                          *
+*******************************************/

#pragma once

#include "glut_program.h"
#include "video_recorder.h"
#include "cps_scene.h"

enum ControlState { kROTATE, kTRANSLATE, kSCALE, kEDIT };

class CPSProgram : public GlutProgram
{
 public:
  CPSProgram();
  virtual ~CPSProgram()
  {
    delete mVideoRecorder;
    delete mScene;
  }

  void Init(int* argc, char* argv[], const char *windowTitle);
  void InitScene(int argc, char *argv[]);

  void DisplayFunc(void);          // Render function.
  void IdleFunc(void);             // Animation/update callback.
  void ReshapeFunc(int w, int h);  // Window resize callback.

  void PassiveMotionFunc(int x, int y);                 // Mouse motion callback.
  void MouseFunc(int button, int state, int x, int y);  // Mouse button callback.
  void MotionFunc(int x, int y);                        // Mouse drag callback.
  void KeyboardFunc(unsigned char key, int x, int y);   // Key pressed.

 private:
  CPSScene* mScene;
  VideoRecorder *mVideoRecorder;

  ControlState mControlState {kROTATE};
};