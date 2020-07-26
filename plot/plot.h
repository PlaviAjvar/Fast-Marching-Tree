#pragma once

#include <vector>
#include <string>
#include <functional>
#include <ostream>
#include "geometry.h"
#include "planning.h"
#include "prm.h"
#include "rrt.h"
#include "fmt.h"

size_t add_graph_edges(
    std::vector <point <double>>& current, 
    std::vector <point <double>>& previous, 
    const std::vector <std::pair<point<double>, point<double>>>& elist
);

size_t add_path_edges(
    std::vector <point <double>>& current, 
    std::vector <point <double>>& previous, 
    const std::vector <std::vector <point <double>>>& paths
);

std::ostream& operator<< (std::ostream& os, std::vector <std::vector <double>> v);

void plot_graph (
    const output& result,
    const std::function<bool(point<double>)>& test_collision,
    const std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&, const std::function<bool(point<double>)>)> add_obstacle_edges,
    const std::vector <std::pair <double, double>> joint_lims,
    bool delay_active = true,
    bool save_image = false
);

std::vector <std::string> get_levels ();

std::string scale_color (
    const std::vector <std::string>& levels,
    const double z,
    const double zlim
);

void plot_3d (
    const output& result,
    const std::function<bool(point<double>)>& test_collision,
    const std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&, const std::function<bool(point<double>)>)> add_obstacle_edges,
    const std::vector <std::pair <double, double>> joint_lims,
    bool delay_active = true,
    bool save_image = false,
    double zlim = 1
);

void disp_snapshot(
    const workspace2d <double>& ws, 
    const point <double> config, 
    const std::string color,
    const std::string name
);

void display_snapshots (
    const workspace2d <double>& ws,
    const std::vector <point <double>>& path
);

