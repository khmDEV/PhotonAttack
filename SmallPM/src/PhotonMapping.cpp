﻿/*********************************************************************************
Copyright (C) 2014 Adrian Jarabo (ajarabo@unizar.es)
Copyright (C) 2014 Diego Gutierrez (diegog@unizar.es)
All rights reserved.

This is an educational Ray Tracer developed for the course 'Informatica Grafica'
(Computer Graphics) tought at Universidad de Zaragoza (Spain). As such, it does not
intend to be fast or general, but just to provide an educational tool for undergraduate
students.

This software is provided as is, and any express or implied warranties are disclaimed.
In no event shall copyright holders be liable for any damage.
**********************************************************************************/
#include "PhotonMapping.h"
#include "World.h"
#include "Intersection.h"
#include "Ray.h"
#include "BSDF.h"
#define _USE_MATH_DEFINES
#include "math.h"



//*********************************************************************
// Compute the photons by tracing the Ray 'r' from the light source
// through the scene, and by storing the intersections with matter
// in the lists 'xx_photons', storing the diffuse (global) and caustic
// photons respectively. For efficiency, both are computed at the same
// time, since computing them separately would result into a lost of
// several samples marked as caustic or diffuse.
// Same goes with the boolean 'direct', that specifies if direct 
// photons (from light to surface) are being stored or not. 
// The initial traced photon has energy defined by the tristimulus
// 'p', that accounts for the emitted power of the light source.
// The function will return true when there are more photons (caustic
// or diffuse) to be shot, and false otherwise.
//---------------------------------------------------------------------
bool PhotonMapping::trace_ray(const Ray& r, const Vector3 &p,
	std::list<Photon> &global_photons, std::list<Photon> &caustic_photons, bool direct)
{

	//Check if max number of shots done...
	if (++m_nb_current_shots > m_max_nb_shots)
	{
		return false;
	}

	// Compute irradiance photon's energy
	Vector3 energy(p);

	Ray photon_ray(r);
	photon_ray.shift();

	bool is_caustic_particle = false;

	//Iterate the path
	while (1)
	{
		// Throw ray and update current_it
		Intersection it;
		world->first_intersection(photon_ray, it);

		if (!it.did_hit())
			break;

		//Check if has hit a delta material...
		if (it.intersected()->material()->is_delta())
		{
			// If delta material, then is caustic...
			// Don't store the photon!
			is_caustic_particle = true;
		}
		else if (photon_ray.get_level() > 0 || direct)
		{
			//If non-delta material, store the photon!
			if (is_caustic_particle)
			{
				//If caustic particle, store in caustics
				if (caustic_photons.size() < m_nb_caustic_photons)
					caustic_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy));
			}
			else
			{
				//If non-caustic particle, store in global
				if (global_photons.size() < m_nb_global_photons)
					global_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy));
			}
			is_caustic_particle = false;
		}

		Real pdf;

		Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
		Real avg_surf_albedo = surf_albedo.avg();

		Real epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);
		while (epsilon2 < 0.)
			epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);

		if (epsilon2 > avg_surf_albedo || photon_ray.get_level() > 5)
			break;

		// Random walk's next step
		// Get sampled direction plus pdf, and update attenuation
		it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf);

		// Shade...
		energy = energy*surf_albedo;
		if (!it.intersected()->material()->is_delta())
			energy *= dot_abs(it.get_normal(), photon_ray.get_direction()) / 3.14159;

		energy = energy / (pdf*avg_surf_albedo);
	}

	if (caustic_photons.size() == m_nb_caustic_photons &&
		global_photons.size() == m_nb_global_photons)
	{
		m_max_nb_shots = m_nb_current_shots - 1;
		return false;
	}

	return true;
}

Vector3 rejectingSampling(){
	Vector3 out;

	do
	{
		double x = (2.0 * static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) - 1.0;
		double y = (2.0 * static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) - 1.0;
		double z = (2.0 * static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) - 1.0;
		out = Vector3(x, y, z);

	} while (out.length() > 1);

	return out;
}

//*********************************************************************
// TODO: Implement the preprocess step of photon mapping,
// where the photons are traced through the scene. To do it,
// you need to follow these steps for each shoot:
//  1 - Sample a world's light source in the scene to create
//		the initial direct photon from the light source.
//	2 - Trace the photon through the scene storing the inter-
//		sections between the photons and matter. You can use
//		the function 'trace_ray' for this purpose.
//	3 - Finally, once all the photons have been shot, you'll
//		need to build the photon maps, that will be used later
//		for rendering. 
//		NOTE: Careful with function
//---------------------------------------------------------------------
void PhotonMapping::preprocess()
{
	bool end = false;

	std::list<Photon> global_photons;
	std::list<Photon> caustic_photons;;

	/*Cargamos los valores de la fuente de luz*/
	Vector3 center = world->light(0).get_position();
	Vector3 Luztotal = world->light(0).get_intensities();

	Vector3 Flux(1, 1, 1);

	/*Se normaliza el flujo de cada foton con respecto a la fuenta de luz*/
	Vector3 FluxNorm = Flux * (Luztotal / (m_nb_caustic_photons + m_nb_global_photons));

	do{

		/*Se genera un rayo aleatorio*/
		Vector3 photonDir = rejectingSampling(); //Vector3 [-1,1]
		Ray photonRay(center, photonDir, 3);

		/*Se lanza el el rayo siguiendo el algorimo normal del trayazor de rayos*/
		end = !trace_ray(photonRay, FluxNorm, global_photons, caustic_photons, false);


	} while (!end);


	/*Se pueblan los KdTrees*/
	for (Photon p : global_photons)
	{
		vector<Real> vec(p.position.data, p.position.data + sizeof(p.position.data) / sizeof(float));
		m_global_map.store(vec, p);
	}
	if (global_photons.size() != 0)
	{
		m_global_map.balance();
	}

	for (Photon p : caustic_photons)
	{
		vector<Real> vec(p.position.data, p.position.data + sizeof(p.position.data) / sizeof(float));
		m_caustics_map.store(vec, p);
	}
	if (caustic_photons.size() != 0)
	{
		m_caustics_map.balance();
	}

}

//*********************************************************************
// TODO: Implement the function that computes the rendering equation 
// using radiance estimation with photon mapping, using the photon
// maps computed as a proprocess. Note that you will need to handle
// both direct and global illumination, together with recursive the 
// recursive evaluation of delta materials. For an optimal implemen-
// tation you should be able to do it iteratively.
// In principle, the class is prepared to perform radiance estimation
// using k-nearest neighbors ('m_nb_photons') to define the bandwidth
// of the kernel.
//---------------------------------------------------------------------


Vector3 PhotonMapping::calculatePhotons(Intersection &it0, bool global, bool caustic) const
{
	Intersection it(it0);

	Real K = 1.0;

	Real position[3] = { it.get_position()[0], it.get_position()[1], it.get_position()[2] };
	const vector<Real> pos(position, position + sizeof(position) / sizeof(Real));
	Vector3 kd = it.intersected()->material()->get_albedo(it);

	Vector3 rada(0.0, 0.0, 0.0);

	if (!m_global_map.is_empty() && global){
		vector<const KDTree<Photon, 3>::Node*> nodes;
		Real max_distance;
		m_global_map.find(pos, m_nb_photons, nodes, max_distance);


		Vector3 flux(0.0, 0.0, 0.0);

		/*Se suma su contribuccion al flujo*/
		for (const KDTree<Photon, 3>::Node* n : nodes)
		{
			Photon pothon = n->data();
			float distance = (pothon.position - it.get_position()).length();
			distance = distance < 0 ? -distance : distance;

			flux += pothon.flux;
		}
		
		Real area = (M_PI*max_distance*max_distance);

		/*Division por el area cubierta*/
		rada = flux.length() == 0 ? Vector3(0.0, 0.0, 0.0) : (flux) / area;

	}
	Vector3 radaCaustic(0.0, 0.0, 0.0);

	if (!m_caustics_map.is_empty() && caustic){
		vector<const KDTree<Photon, 3>::Node*> nodes;
		Real max_distance;

		Vector3 flux(0.0, 0.0, 0.0);

		m_caustics_map.find(pos, m_nb_photons, nodes, max_distance);

		/*Se suma su contribuccion al flujo*/
		for (const KDTree<Photon, 3>::Node* n : nodes)
		{
			Photon pothon = n->data();
			float distance = (pothon.position - it.get_position()).length();
			distance = distance < 0 ? -distance : distance;

			flux += pothon.flux;
		}

		Real area = (M_PI*max_distance*max_distance);

		/*Division por el area cubierta*/
		radaCaustic = flux.length() == 0 ? Vector3(0.0, 0.0, 0.0) : (flux) / area;

	}

	Vector3 res = rada + radaCaustic;

	return res;
}

Vector3 PhotonMapping::calculateDirect(Intersection &it0) const
{
	Intersection it(it0);
	Vector3 res(0);

	Real pdf;

	for (int o = 0; o < 5; o++)
	{
		Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
		Ray photon_ray = it.get_ray();


		if (it.intersected()->material()->is_delta())
		{
			it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf);
			photon_ray.shift();
			world->first_intersection(photon_ray, it);
		}
		else{
			for (int i = 0; i < world->nb_lights(); i++)
			{
				if (world->light(i).is_visible(it.get_position())){
					Vector3 lightDirection = normalize(world->light(i).get_position() - it.get_position());
					Real angle = dot(it.get_normal(), lightDirection);
					angle = angle < 0 ? 0 : angle;
					Vector3 light = world->light(i).get_incoming_light(it.get_position())*(angle);
					res = light*surf_albedo;
				}

			}
			break;
		}

	}

	return res;
}

Vector3 PhotonMapping::shade(Intersection &it0)const
{

	Vector3 L(0);
	Intersection it(it0);



	//**********************************************************************
	// The following piece of code is included here for two reasons: first
	// it works as a 'hello world' code to check that everthing compiles 
	// just fine, and second, to illustrate some of the functions that you 
	// will need when doing the work. Goes without saying: remove the 
	// pieces of code that you won't be using.
	//

	unsigned int debug_mode = 0;
	Vector3 directa;
	Vector3 global;
	Vector3 caustica;

	switch (debug_mode)
	{
	case 0:
		// ----------------------------------------------------------------
		// Photons
		directa = calculateDirect(it);

		global = calculatePhotons(it, true, false);

		caustica = calculatePhotons(it, false, true);


		L = directa + global + caustica;

		//L = L.length() == 0 ? L : L.normalize();

		break;
	case 7:
		// ----------------------------------------------------------------
		// Photons
		L = calculatePhotons(it, true, true);


		break;
	case 8:
		// ----------------------------------------------------------------
		// Photons
		L = calculateDirect(it);

		break;
	case 1:
		// ----------------------------------------------------------------
		// Display Albedo Only
		L = it.intersected()->material()->get_albedo(it);
		break;
	case 2:
		// ----------------------------------------------------------------
		// Display Normal Buffer
		L = it.get_normal();
		break;
	case 3:
		// ----------------------------------------------------------------
		// Display whether the material is specular (or refractive) 
		L = Vector3(it.intersected()->material()->is_delta());
		break;

	case 4:
		// ----------------------------------------------------------------
		// Display incoming illumination from light(0)
		L = world->light(0).get_incoming_light(it.get_position());
		break;

	case 5:
		// ----------------------------------------------------------------
		// Display incoming direction from light(0)
		L = world->light(0).get_incoming_direction(it.get_position());
		break;

	case 6:
		// ----------------------------------------------------------------
		// Check Visibility from light(0)
		if (world->light(0).is_visible(it.get_position()))
			L = Vector3(1.);
		break;
	}
	// End of exampled code
	//**********************************************************************

	return L;
}