#pragma once

#include <vector>
#include <functional>
#include <string>
#include <functional>
#include <ostream>
#include "geometry.h"
#include "planning.h"
#include "prm.h"
#include "rrt.h"
#include "fmt.h"
#include <map>
#include "matplotlibcpp.h"
#include <random>
namespace plt = matplotlibcpp;

class plot_object {
private:
    std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&)> add_obstacle_edges;
    std::vector <std::pair <double, double>> joint_limits;
    std::function<bool(point<double>)> test_collision;
    workspace3d <double> ws;
    std::string mode;
    unsigned int resolution;

    // rotate around z-axis
    void rotate (
        double& x,
        double& y,
        double phi
    ) const {
        double nx = x * cos(phi) + y * sin(phi);
        double ny = -x * sin(phi) + y * cos(phi);
        x = nx;
        y = ny;
    }

public:
    plot_object () {
        mode = "normal"; 
    }

    plot_object (
        const plot_object& po
    ) : add_obstacle_edges(po.add_obstacle_edges), joint_limits(po.joint_limits), test_collision(po.test_collision), 
        mode(po.mode), resolution(po.resolution) {}

    plot_object (
        const std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&)> aoe
    ) : add_obstacle_edges(aoe) {
        mode = "normal";
    }

    plot_object (
        const std::vector <std::pair <double, double>>& lims,
        const std::function<bool(point<double>)> tc,
        const unsigned int res = 200
    ) : joint_limits(lims), test_collision(tc), resolution(res) {
        mode = "scatter";
    }

    plot_object (
        const workspace3d <double>& _ws
    ) : ws(_ws) {
        mode = "normal3d";
    }

    void plot (
        const std::string& flag
    ) const {

        // normal plotting (plot edges of polygon)
        if (mode == "normal") {
            std::vector <point <double>> cur, prev;
            size_t cutoff = add_obstacle_edges(cur, prev);

            for (size_t i = 0; i < cutoff; ++i) {
                std::vector <double> xs{prev[i][0], cur[i][0]};
                std::vector <double> ys{prev[i][1], cur[i][1]};
                plt::plot(xs, ys, flag);
            }

            return;
        }

        // scatter plot to construct configuration space obstacles
        if (mode == "scatter") { 
            if (joint_limits.size() != 2) {
                throw std::length_error("Manipulator has to be 2d for this to work");
            }
            double dx = (joint_limits[0].second - joint_limits[0].first) / resolution;
            double dy = (joint_limits[1].second - joint_limits[1].first) / resolution;
            std::vector <double> xs, ys;

            for (size_t yi = 1; yi < resolution; ++yi) {
                for (size_t xi = 1; xi < resolution; ++xi) {
                    double x = joint_limits[0].first + xi * dx;
                    double y = joint_limits[1].first + yi * dy;
                    point <double> cfg = point2d<double>(x,y);

                    // test configuration for collision
                    if (test_collision(cfg)) {
                        xs.push_back(x);
                        ys.push_back(y);
                    }
                }
            }

            std::map <std::string, std::string> flm;
            flm["color"] = flag[0];
            flm["linestyle"] = flag[1];
            plt::scatter(xs, ys, 1.0, flm);
            return;
        }

        // never should reach here
        throw std::domain_error("Mode has to be normal or scatter");
    }


    void plot3d (
        std::vector <double>& xs,
        std::vector <double>& ys,
        std::vector <double>& zs,
        const double eps = 0.1
    ) const {
        // 3d normal plotting mode
        if (mode == "normal3d") {
            for (const box <double>& abox : ws.get_obstacles()) {
                std::vector <std::pair <point <double>, point <double>>> edges = abox.get_edges();
                
                for (const auto& edge : edges) {
                    point3d <double> fi = edge.first;
                    point3d <double> se = edge.second;
                    xs.push_back(fi.getx());
                    xs.push_back(se.getx());
                    ys.push_back(fi.gety());
                    ys.push_back(se.gety());
                    zs.push_back(fi.getz());
                    zs.push_back(se.getz());
                }

                xs.push_back(NAN);
                ys.push_back(NAN);
                zs.push_back(NAN);
            }
            return;
        }

        // scatter plotting
        // have to
        if (mode == "scatter") { 
            if (joint_limits.size() != 3) {
                throw std::length_error("Manipulator has to be 3d for this to work");
            }
            double dx = (joint_limits[0].second - joint_limits[0].first) / resolution;
            double dy = (joint_limits[1].second - joint_limits[1].first) / resolution;
            double dz = (joint_limits[2].second - joint_limits[2].first) / resolution;

            std::default_random_engine gen;
            std::uniform_real_distribution <double> dist(0, 2*M_PI);

            for (size_t zi = 1; zi < resolution; ++zi) {
                for (size_t yi = 1; yi < resolution; ++yi) {
                    for (size_t xi = 1; xi < resolution; ++xi) {
                        double x = joint_limits[0].first + xi * dx;
                        double y = joint_limits[1].first + yi * dy;
                        double z = joint_limits[2].first + zi * dz;
                        point <double> cfg = point3d<double>(x,y,z);

                        // test configuration for collision
                        if (test_collision(cfg)) {
                            xs.push_back(x);
                            ys.push_back(y);
                            zs.push_back(z);

                            // construct box in obstacle space
                            for (int xd = 0; xd <= 1; ++xd) {
                                for (int yd = 0; yd <= 1; ++yd) {
                                    for (int zd = 0; zd <= 1; ++zd) {
                                        double nx = eps * xd;
                                        double ny = eps * yd;
                                        double nz = eps * zd;

                                        xs.push_back(x + nx);
                                        ys.push_back(y + ny);
                                        zs.push_back(z + nz);
                                    }
                                }
                            }

                            xs.push_back(NAN);
                            ys.push_back(NAN);
                            zs.push_back(NAN);
                        }
                    }
                }
            }
            return;
        }


        // never should reach here
        throw std::domain_error("Mode has to be normal3d or scatter");
    }
};

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
    const plot_object& pltobj,
    const std::vector <std::pair <double, double>> joint_lims,
    bool delay_active = true,
    bool save_image = false
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

void plot3d (
    const output& result,
    const plot_object& pltobj,
    const std::vector <std::pair <double, double>> joint_lims,
    bool show_path = false,
    bool save_image = false
);

void disp_snapshot3d (
    const workspace3d <double>& ws, 
    const point <double> config, 
    std::vector <double>& xs,
    std::vector <double>& ys,
    std::vector <double>& zs
);

void display_snapshots3d (
    const workspace3d <double>& ws,
    const std::vector <point <double>>& path
);

void plot_function (
    const std::vector <double>& xs,
    const std::vector <double>& ys,
    const std::string xlabel,
    const std::string ylabel,
    const std::string title
);

void obtain_limits (
    const workspace2d <double>& ws,
    const std::vector <point <double>>& path,
    double& xlow,
    double& xhigh,
    double& ylow,
    double& yhigh,
    const double offset = 0.2
);