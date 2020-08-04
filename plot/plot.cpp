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
    const plot_object& pltobj,
    const std::vector <std::pair <double, double>> joint_lims,
    bool delay_active,
    bool save_image
) {
    plt::backend("TkAgg");
    std::vector <point <double>> current, previous;

    // add obstacle edges
    pltobj.plot("g-");

    // add graph edges
    size_t cutoffB = add_graph_edges(current, previous, result.get_edgelist());

    // add path edges
    size_t cutoffC = add_path_edges(current, previous, result.get_paths());

    plt::xlim(joint_lims[0].first, joint_lims[0].second);
    plt::ylim(joint_lims[1].first, joint_lims[1].second);

    for (size_t i = 0; i < cutoffC; ++i) {
        std::string lineflag;

        if (i < cutoffB) lineflag = "b-";
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


void obtain_limits (
    const std::vector <polygon <double>>& obs,
    double& xlow,
    double& xhigh,
    double& ylow,
    double& yhigh,
    const double offset
) {

    // find leftmost/rightmost x-coordinate and bottommost/topmost y-coordinate
    for (const auto& polygon : obs) {
        for (const auto& edge : polygon.get_edges()) {
            point2d <double> prev = edge.first;
            point2d <double> cur = edge.second;
            xlow = std::min(xlow, cur.getx());
            xhigh = std::max(xhigh, cur.getx());
            ylow = std::min(ylow, cur.gety());
            yhigh = std::max(yhigh, cur.gety());
        }
    }

    // add offset for visualization
    xlow -= offset;
    xhigh += offset;
    ylow -= offset;
    yhigh += offset;

    // adjust equal bounds
    double dif = std::max(xhigh - xlow, yhigh - ylow);
    const double epsilon = 1e-4;

    if (xhigh - xlow < dif - epsilon) {
        double half = (dif - (xhigh - xlow)) / 2;
        xlow -= half;
        xhigh += half;
    }

    if (yhigh - ylow < dif - epsilon) {
        double half = (dif - (yhigh - ylow)) / 2;
        ylow -= half;
        yhigh += half;
    }
}


// display snapshots of 2D planar arm in its movement from start to goal
void display_snapshots (
    const workspace2d <double>& ws,
    const std::vector <point <double>>& path
) {
    plt::backend("TkAgg");
    std::vector <std::string> names{"start","2nd","3rd","4th","5th","goal"};
    std::vector <std::string> colors{"y", "m", "c", "r", "b", "k"};
    size_t num_snapshots = colors.size();

    // obtain xy-limits
    double xlow, xhigh, ylow, yhigh;
    obtain_limits(ws.get_obstacles(), xlow, xhigh, ylow, yhigh);

    plt::figure_size(500, 500);
    plt::xlim(xlow, xhigh);
    plt::ylim(ylow, yhigh);

    for (const auto& polygon : ws.get_obstacles()) {
        for (const auto& edge : polygon.get_edges()) {
            point2d <double> prev = edge.first;
            point2d <double> cur = edge.second;
            plt::plot(std::vector <double>{prev.getx(), cur.getx()}, std::vector <double>{prev.gety(), cur.gety()}, "g-");
        }
    }

    if (path.empty()) {
        std::cout << "Path empty" << std::endl;
        return;
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
    const plot_object& pltobj,
    const std::vector <std::pair <double, double>> joint_lims,
    bool show_path,
    bool save_image
) {
    plt::backend("TkAgg");
    std::vector <point <double>> current, previous;
    std::vector <double> xs, ys, zs;

    // plot the obstacles
    pltobj.plot3d(xs, ys, zs);

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

    if (path.empty()) {
        std::cout << "Path empty" << std::endl;
        return;
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

void plot_function (
    const std::vector <double>& xs,
    const std::vector <double>& ys,
    const std::string xlabel,
    const std::string ylabel,
    const std::string title
) {
    plt::backend("TkAgg");
    plt::plot(xs, ys);
    plt::title(title);
    plt::xlabel(xlabel);
    plt::ylabel(ylabel);
    plt::show();
}