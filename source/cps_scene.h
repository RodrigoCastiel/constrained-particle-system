/******************************************+
*                                          *
*  Author: Rodrigo Castiel                 *
*  Email: rcrs2@cin.ufpe.br                *
*         castielr@usc.edu                 *
*                                          *
+*******************************************/

#pragma once

#include "scene.h"
#include "particle_system.h"

class CPSScene : public Scene
{
public:
  CPSScene() { }

  void Init(BasicPipelineProgram* pipelineProgram, GLuint programHandle);
  void Clean();

  void Render();
  void Animate();

  void OnMouseLeftClick(int x, int y, int w, int h, int state);
  void OnMouseRightClick(int x, int y, int w, int h);

  void OnMouseLeftDrag(int x, int y, int w, int h);

  void TogglePauseSimulation() { mParticleSystem->TogglePauseSimulation(); }

  virtual ~CPSScene() { Scene::Clean(); }

public:
  ParticleSystem* mParticleSystem;

};
