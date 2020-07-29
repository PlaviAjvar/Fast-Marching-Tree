#include "plot.h"
#include "matplotlibcpp.h"
#include "geometry.h"
#include <sstream>
namespace plt = matplotlibcpp;

size_t add_graph_edges(
    std::vector <point <double>>& current, 
    std::vector <point <double>>& previous, 
    const std::vector <std::pair<point<double>, point<double>>>& elist
) {

    for (const auto& vertex : elist) {
        point2d <double> prev_point(vertex.first);
        point2d <double> cur_point(vertex.second);

        previous.push_back(prev_point);
        current.push_back(cur_point);
    }

    return current.size();
}


size_t add_path_edges(
    std::vector <point <double>>& current, 
    std::vector <point <double>>& previous, 
    const std::vector <std::vector <point <double>>>& paths
) {

    for (auto&& path : paths) {
        for (size_t i = 1; i < path.size(); ++i) {
            previous.push_back(path[i-1]);
            current.push_back(path[i]);
        }
    }

    return current.size();
}

std::ostream& operator<< (std::ostream& os, std::vector <std::vector <double>> v) {
    for (auto row : v) {
        for (auto el : v) {
            os << el << " ";
        }
        os << std::endl;
    }
    return os;
}

// plot 2D graph
void plot_graph (
    const output& result,
    const std::function<bool(point<double>)>& test_collision,
    const std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&, const std::function<bool(point<double>)>)> add_obstacle_edges,
    const std::vector <std::pair <double, double>> joint_lims,
    bool delay_active,
    bool save_image
) {
    plt::backend("TkAgg");
    std::vector <point <double>> current, previous;

    // add obstacle edges
    size_t cutoffA = add_obstacle_edges(current, previous, test_collision);

    // add graph edges
    size_t cutoffB = add_graph_edges(current, previous, result.get_edgelist());

    // add path edges
    size_t cutoffC = add_path_edges(current, previous, result.get_paths());

    plt::xlim(joint_lims[0].first, joint_lims[0].second);
    plt::ylim(joint_lims[1].first, joint_lims[1].second);

    for (size_t i = 0; i < cutoffC; ++i) {
        std::string lineflag;

        if (i < cutoffA) lineflag = "g-";
        else if (i < cutoffB) lineflag = "b-";
        else lineflag = "r-";

        point2d <double> prev(previous[i]);
        point2d <double> cur(current[i]);
        plt::plot(std::vector <double>{prev.getx(), cur.getx()}, std::vector <double>{prev.gety(), cur.gety()}, lineflag); 
        
        if (delay_active) {
            plt::pause(0.001);
        }
    }

    if (!save_image) {
        plt::show();
    }
    else {
        plt::save("./output.png");
    }
}

// get as many as possible discrete levels of green-blue
std::vector <std::string> get_levels () {
    std::vector <std::string> levels;

    for (unsigned int lev = 1; lev < 256; ++lev) {
        std::stringstream ss;
        ss << std::hex << lev;
        std::string res = ss.str();

        for(auto &ch : res) {
            ch = toupper(ch);
        }
        if (res.size() == 1) res = "0" + res;

        levels.push_back("#00" + res + "FF");
    }

    return levels;
}

// get discrete level of blue from z coordinate
std::string scale_color (
    const std::vector <std::string>& levels,
    const double z,
    const double zlim
) {
    size_t levelcnt = levels.size();
    double quant = zlim / levelcnt;
    size_t lev = z / quant;
    return levels[lev];
}

// plot 3D graph
void plot_bird (
    const output& result,
    const std::function<bool(point<double>)>& test_collision,
    const std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&, const std::function<bool(point<double>)>)> add_obstacle_edges,
    const std::vector <std::pair <double, double>> joint_lims,
    bool delay_active,
    bool save_image,
    double zlim
) {
    plt::backend("TkAgg");
    std::vector <point <double>> current, previous;

    // add obstacle edges
    size_t cutoffA = add_obstacle_edges(current, previous, test_collision);

    // add graph edges
    size_t cutoffB = add_graph_edges(current, previous, result.get_edgelist());

    // add path edges
    size_t cutoffC = add_path_edges(current, previous, result.get_paths());

    plt::xlim(joint_lims[0].first, joint_lims[0].second);
    plt::ylim(joint_lims[1].first, joint_lims[1].second);

    std::vector <std::string> levels(get_levels());

    for (size_t i = 0; i < cutoffC; ++i) {
        std::string lineflag;
        point3d <double> prev(previous[i]);
        point3d <double> cur(current[i]);

        if (i < cutoffA) {
            lineflag = "g-";
        }
        else if (i < cutoffB) {
            // if blue scale it
            lineflag = scale_color(levels, cur.getz(), zlim);
        }
        else {
            lineflag = "r-";
        }

        plt::plot(std::vector <double>{prev.getx(), cur.getx()}, std::vector <double>{prev.gety(), cur.gety()}, lineflag); 
        
        if (delay_active) {
            plt::pause(0.001);
        }
    }

    if (!save_image) {
        plt::show();
    }
    else {
        plt::save("./output.png");
    }
}

// display single snapshot
void disp_snapshot(
    const workspace2d <double>& ws, 
    const point <double> config, 
    const std::string color,
    const std::string name
) {
    arm <double>& robot = ws.get_robot();
    auto edges = robot.dir_kine(config);
    std::vector <double> xs;
    std::vector <double> ys;

    for (const auto& space_edge : edges) {
        std::pair <point2d <double>, point2d <double>> edge(space_edge);
        xs.push_back(edge.first.getx());
        xs.push_back(edge.second.getx());
        ys.push_back(edge.first.gety());
        ys.push_back(edge.second.gety());
    }

    plt::named_plot(name, xs, ys, color + "-");
}


// display snapshots of 2D planar arm in its movement from start to goal
void display_snapshots (
    const workspace2d <double>& ws,
    const std::vector <point <double>>& path
) {
    plt::backend("TkAgg");
    std::vector <std::string> names{"start","2nd","3rd","4th","5th","goal"};
    std::vector <std::string> colors{"y", "m", "c", "r", "g", "b"};
    size_t num_snapshots = colors.size();

    for (const auto& polygon : ws.get_obstacles()) {
        for (const auto& edge : polygon.get_edges()) {
            point2d <double> prev = edge.first;
            point2d <double> cur = edge.second;
            plt::plot(std::vector <double>{prev.getx(), cur.getx()}, std::vector <double>{prev.gety(), cur.gety()}, "g-");
        }
    }

    disp_snapshot(ws, path[0], colors[0], names[0]);

    for (size_t i = 1; i < num_snapshots - 1; ++i) {
        double jump = double(path.size()) / num_snapshots;
        size_t path_idx = i * jump;
        disp_snapshot(ws, path[path_idx], colors[i], names[i]);
    }

    disp_snapshot(ws, path.back(), colors.back(), names.back());
    plt::legend();
    plt::show();
}

std::ostream& operator<< (std::ostream& os, std::vector <double> v) {
    for (auto e : v) os << e << " ";
    return os;
}

// plot in 3d
void plot3d (
    const output& result,
    const workspace3d <double>& ws,
    const std::vector <std::pair <double, double>> joint_lims,
    bool show_path,
    bool save_image,
    bool show_obstacles
) {
    plt::backend("TkAgg");
    std::vector <point <double>> current, previous;
    std::vector <double> xs, ys, zs;

    // plot the obstacles
    if (show_obstacles) {
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
    }

    // add graph edges or path edges
    size_t cutoff;
    
    if (!show_path) {
        cutoff = add_graph_edges(current, previous, result.get_edgelist());
    } 
    else {
        cutoff = add_path_edges(current, previous, result.get_paths());
    }

    for (size_t i = 0; i < cutoff; ++i) {
        point3d <double> fi(previous[i]);
        point3d <double> se(current[i]);

        xs.push_back(fi.getx());
        xs.push_back(se.getx());
        ys.push_back(fi.gety());
        ys.push_back(se.gety());
        zs.push_back(fi.getz());
        zs.push_back(se.getz());
        xs.push_back(NAN);
        ys.push_back(NAN);
        zs.push_back(NAN);
    }

    plt::plot3(xs, ys, zs);

    if (!save_image) {
        plt::show();
    }
    else {
        plt::save("./output.png");
    }
}


// display single snapshot of 3d arm
void disp_snapshot3d (
    const workspace3d <double>& ws, 
    const point <double> config, 
    std::vector <double>& xs,
    std::vector <double>& ys,
    std::vector <double>& zs
) {
    arm <double>& robot = ws.get_robot();
    auto edges = robot.dir_kine(config);

    for (const auto& space_edge : edges) {
        std::pair <point3d <double>, point3d <double>> edge = linkto3d(space_edge);
        std::cout << edge.first << " " << edge.second << std::endl;
        xs.push_back(edge.first.getx());
        xs.push_back(edge.second.getx());
        ys.push_back(edge.first.gety());
        ys.push_back(edge.second.gety());
        zs.push_back(edge.first.getz());
        zs.push_back(edge.second.getz());
    }

    xs.push_back(NAN);
    ys.push_back(NAN);
    zs.push_back(NAN);
}


// dsplay snapshots of 3d arm
void display_snapshots3d (
    const workspace3d <double>& ws,
    const std::vector <point <double>>& path
) {
    plt::backend("TkAgg");
    size_t num_snapshots = 10;
    std::vector <double> xs, ys, zs;

    // plot the obstacles
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
    
    disp_snapshot3d(ws, path[0], xs, ys, zs);

    for (size_t i = 1; i < num_snapshots - 1; ++i) {
        double jump = double(path.size()) / num_snapshots;
        size_t path_idx = i * jump;
        disp_snapshot3d(ws, path[path_idx], xs, ys, zs);
    }

    disp_snapshot3d(ws, path.back(), xs, ys, zs);

    plt::plot3(xs, ys, zs);
    plt::show();
}
