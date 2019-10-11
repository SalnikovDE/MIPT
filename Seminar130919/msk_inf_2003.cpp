#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <exception>
#include <cstdlib>
#include <string>

/*
 Карта мира в компьютерной игре “Цивилизация” версии 1 представляет собой прямоугольник, разбитый на квадратики.
 Каждый квадратик может иметь один из нескольких возможных рельефов, для простоты ограничимся тремя видами
 рельефов - поле, лес и вода. Поселенец перемещается по карте, при этом на перемещение в клетку, занятую полем,
 необходима одна единица времени, на перемещение в лес - две единицы времени, а перемещаться в клетку с водой нельзя.

У вас есть один поселенец, вы определили место, где нужно построить город, чтобы как можно скорее завладеть всем миром.
 Найдите маршрут переселенца, приводящий его в место строительства города, требующий минимального времени. Н
 а каждом ходе переселенец может перемещаться в клетку, имеющую общую сторону с той клеткой, где он сейчас находится.
*/


class Graph {
public:
    typedef size_t Vertex;

    Graph(size_t _vertex_count, bool _is_directed) {
        is_directed = _is_directed;
        vertex_count = _vertex_count;
        edges_count = 0;
    }

    bool directed() {
        return is_directed;
    }

    virtual size_t get_vertex_deg(const Vertex &v) const = 0;

    size_t get_vertex_count() const {
        return vertex_count;
    }

    size_t get_edges_count() const {
        return edges_count;
    }

    virtual void add_edge(const Vertex &start, const Vertex &finish) {
        edges_count++;
    }

    virtual std::vector<Vertex> get_all_neighbours(const Vertex &v) const = 0;

protected:
    size_t vertex_count, edges_count;
    bool is_directed;
};

class AdjListGraph : public Graph {
public:
    AdjListGraph(size_t _vertex_count, bool _is_directed) : Graph(_vertex_count, _is_directed) {
        Vertexes.resize(vertex_count);
    }

    AdjListGraph(bool _is_directed, const std::vector<std::vector<int>> &matrix) : Graph(
            matrix.size(), _is_directed) {
        Vertexes.resize(vertex_count);
        read_matrix(matrix);
    }

    virtual void add_edge(const Vertex &start, const Vertex &finish) override {
        Graph::add_edge(start, finish);
        Vertexes[start].push_back(finish);
        if (!is_directed) {
            Vertexes[finish].push_back(start);
        }
    }

    virtual size_t get_vertex_deg(const Vertex &v) const override {
        return Vertexes[v].size();
    }


    virtual std::vector<Vertex> get_all_neighbours(const Vertex &v) const override {
        return Vertexes[v];
    }

private:
    std::vector<std::vector<Vertex>> Vertexes;

    void read_matrix(const std::vector<std::vector<int>> &matrix) {
        if (matrix.size() != get_vertex_count()) {
            throw std::out_of_range("wrong matrix size");
        }
        for (Vertex i = 0; i < matrix.size(); i++) {
            if (matrix[i].size() != get_vertex_count()) {
                throw std::out_of_range("wrong matrix size");
            }
            for (Vertex j = 0; j < matrix[i].size(); j++) {
                if (matrix[i][j] != 0) {
                    if (j >= i || is_directed) {
                        add_edge(i, j);
                    }
                }
            }
        }
    };
};

class WeightedAgjListGraph : public AdjListGraph {
private:
    std::vector<int> weight;
public:

};
namespace GraphAlgorithms {
    enum vertexColours {
        WHITE, GRAY, BLACK
    };
    typedef std::vector<Graph::Vertex> Component;

    void dfs_graph_top_sort(Graph::Vertex ver, const Graph &graph, std::vector<vertexColours> &colours,
                            std::vector<Graph::Vertex> &result) {
        colours[ver] = GRAY;
        for (const Graph::Vertex &v : graph.get_all_neighbours(ver)) {
            if (colours[v] == WHITE) {
                dfs_graph_top_sort(v, graph, colours, result);
            }
        }
        colours[ver] = BLACK;
        result.push_back(ver);

    }

    void dfs_find_components(const Graph::Vertex &v, const Graph &graph,
                             std::vector<vertexColours> &colours, Component &component) {
        colours[v] = GRAY;
        component.push_back(v);
        for (const Graph::Vertex &i : graph.get_all_neighbours(v)) {
            if (colours[i] == WHITE) {
                dfs_find_components(i, graph, colours, component);
            }
        }
        colours[v] = BLACK;
    }

    bool dfs_check_cicle(const Graph::Vertex &v, const Graph &graph,
                         std::vector<vertexColours> &colours) {
        colours[v] = GRAY;
        for (const Graph::Vertex &i : graph.get_all_neighbours(v)) {
            if (colours[i] == WHITE) {
                if (dfs_check_cicle(i, graph, colours)) return true;
            } else if (colours[i] == GRAY) {
                colours[v] = BLACK;
                return true;
            }
        }
        colours[v] = BLACK;
        return false;
    }

    std::vector<Component> find_components(const Graph &graph) {
        std::vector<Component> result;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                result.emplace_back(std::vector<Graph::Vertex>());
                dfs_find_components(i, graph, colours, result.back());
            }
        }

        return result;
    }

    bool is_cicled(const Graph &graph) {
        bool result = 0;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE && !result) {
                result = dfs_check_cicle(i, graph, colours);
            }
        }
        return result;
    }

    std::vector<Graph::Vertex>
    find_shortest_way_bfs(const Graph &graph, const Graph::Vertex &start, const Graph::Vertex &finish) {
        const int NOT_SET = -1;
        std::queue<Graph::Vertex> to_visit;
        std::vector<Graph::Vertex> ancestor(graph.get_vertex_count());
        std::vector<int> distance(graph.get_vertex_count(), NOT_SET);
        to_visit.push(start);
        distance[start] = 0;

        while (!to_visit.empty()) {
            auto cur = to_visit.front();

            to_visit.pop();
            for (const Graph::Vertex &n : graph.get_all_neighbours(cur)) {
                if (distance[n] == NOT_SET) {
                    distance[n] = distance[cur] + 1;
                    ancestor[n] = cur;
                    to_visit.push(n);
                }
            }
        }
        std::vector<Graph::Vertex> result;
        auto cur = finish;
        for (int i = 0; i < distance[finish] + 1; i++) {
            result.push_back(cur);
            cur = ancestor[cur];
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    std::vector<Graph::Vertex> graph_top_sort(const Graph &graph) {
        std::vector<Graph::Vertex> result;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                dfs_graph_top_sort(i, graph, colours, result);
            }
        }

        return result;
    }


};

void
add(Graph &g, int i_start, int j_start, int i_finish, int j_finish, const std::vector<std::vector<char>> &world_map,
    int &fake) {

    if (world_map[i_finish][j_finish] == '#') {
        return;
    } else if (world_map[i_finish][j_finish] == '.') {
        g.add_edge(i_start * world_map[0].size() + j_start, i_finish * world_map[0].size() + j_finish);
    } else if (world_map[i_finish][j_finish] == 'W') {
        g.add_edge(i_start * world_map[0].size() + j_start, fake);
        g.add_edge(fake, i_finish * world_map[0].size() + j_finish);
        fake++;
    }
}

char get_direction(int a, int b, int cnt_columns) {
    if (b - a == cnt_columns) {
        return 'S';
    }
    if (b - a == -cnt_columns) {
        return 'N';
    }
    if (b - a == 1) {
        return 'E';
    }
    if (b - a == -1) {
        return 'W';
    }
}

int get_graph_vertex(int i, int j, int cnt_columns) {
    return j + i * cnt_columns;
}

int main() {
    int cnt_raws, cnt_columns;
    int i_start, j_start;
    int i_finish, j_finish;

    std::cin >> cnt_raws >> cnt_columns;
    std::cin >> i_start >> j_start;
    std::cin >> i_finish >> j_finish;

    i_start--;
    i_finish--;
    j_start--;
    j_finish--;

    int max_real_vertex = cnt_raws * cnt_columns - 1;
    int woods_count = 0;
    std::vector<std::vector<char>> world_map(cnt_raws, std::vector<char>(cnt_columns));

    for (auto &s : world_map) {
        for (auto &p : s) {
            std::cin >> p;
            if (p == 'W') {
                woods_count++;
            }
        }
    }

    AdjListGraph g(max_real_vertex + 4 * woods_count + 1, true);

    int fake = max_real_vertex + 1;
    for (int i = 0; i < cnt_raws; i++) {
        for (int j = 0; j < cnt_columns; j++) {
            if (i + 1 < cnt_raws) {
                add(g, i, j, i + 1, j, world_map, fake);
            }
            if (i - 1 >= 0) {
                add(g, i, j, i - 1, j, world_map, fake);
            }
            if (j + 1 < cnt_columns) {
                add(g, i, j, i, j + 1, world_map, fake);
            }
            if (j - 1 >= 0) {
                add(g, i, j, i, j - 1, world_map, fake);
            }
        }
    }
    auto result = GraphAlgorithms::find_shortest_way_bfs(g, get_graph_vertex(i_start, j_start, cnt_columns),
                                                         get_graph_vertex(i_finish, j_finish, cnt_columns));
    if (result.size() == 0) {
        std::cout << -1;
    } else {
        std::cout << result.size() - 1 << std::endl;
        for (int i = 0; i < result.size() - 1; i++) {
            if (result[i + 1] > max_real_vertex) {
                std::cout << get_direction(result[i], result[i + 2], cnt_columns);
                i++;
            } else {
                std::cout << get_direction(result[i], result[i + 1], cnt_columns);
            }
        }
    }
}