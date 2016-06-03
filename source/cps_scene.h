/******************************************+
*                                          *
*  Author: Rodrigo Castiel                 *
*  Email: rcrs2@cin.ufpe.br                *
*         castielr@usc.edu                 *
*                                          *
+*******************************************/

#pragma once

#include "scene.h"

class CPSScene : public Scene
{
public:
  CPSScene() { }

  void Init(BasicPipelineProgram* pipelineProgram, GLuint programHandle);
  void Clean();

  void Render();
  void Animate();

  void OnMouseLeftClick(int x, int y, int w, int h);
  void OnMouseRightClick(int x, int y, int w, int h);

  virtual ~CPSScene() { Scene::Clean(); }
};
