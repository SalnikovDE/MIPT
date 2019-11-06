#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <exception>
#include <cstdlib>
#include <string>
#include <set>
#include <unordered_map>

/*
 * Во время контрольной работы профессор Флойд заметил, что некоторые студенты обмениваются записками. Сначала он хотел поставить им всем двойки, но в тот день профессор был добрым, '
     а потому решил разделить студентов на две группы: списывающих и дающих списывать, и поставить двойки только первым.
     У профессора записаны все пары студентов, обменявшихся записками. Требуется определить, сможет ли он разделить студентов на две группы так, 
     чтобы любой обмен записками осуществлялся от студента одной группы студенту другой группы.
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

class WeightedGraph {
public:

};

class WeightedAgjListGraph : public AdjListGraph, public WeightedGraph {
private:
    std::vector<int> weight;
public:

};

namespace GraphAlgorithms {
    enum vertexColours {
        WHITE, GRAY, BLACK
    };
    typedef std::vector<Graph::Vertex> Component;
    typedef std::pair<Graph::Vertex, Graph::Vertex> Edge;


    void dfs_graph_top_sort(Graph::Vertex ver, const Graph &graph, std::vector<vertexColours> &colours,
                            std::vector<Graph::Vertex> &top_sort_result) {
        colours[ver] = GRAY;
        for (const Graph::Vertex &v : graph.get_all_neighbours(ver)) {
            if (colours[v] == WHITE) {
                dfs_graph_top_sort(v, graph, colours, top_sort_result);
            }
        }
        colours[ver] = BLACK;
        top_sort_result.push_back(ver);

    }

    bool dfs_bipartite(const Graph::Vertex &v, const Graph &graph, std::vector<vertexColours> &colours,
                       std::vector<int> &part) {
        colours[v] = GRAY;
        for (const Graph::Vertex &i : graph.get_all_neighbours(v)) {
            if (colours[i] == GRAY) {
                if (part[i] == part[v]) {
                    return false;
                }
            } else if (colours[i] == WHITE) {
                part[i] = -part[v];
                if(!dfs_bipartite(i, graph, colours, part)) {
                    return false;
                }
            }
        }
        colours[v] = BLACK;
        return true;
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


    Graph::Vertex dfs_find_cicle(const Graph::Vertex &v, const Graph &graph,
                                 std::vector<vertexColours> &colours, std::vector<Graph::Vertex> &ancestors) {
        colours[v] = GRAY;
        Graph::Vertex last_vertex;
        const Graph::Vertex NOT_FOUND = graph.get_vertex_count();

        for (const Graph::Vertex &i : graph.get_all_neighbours(v)) {
            if (colours[i] == GRAY) {
                colours[v] = BLACK;
                ancestors[i] = v;
                return i;
            } else if (colours[i] == WHITE) {
                ancestors[i] = v;
                last_vertex = dfs_find_cicle(i, graph, colours, ancestors);
                if (last_vertex != NOT_FOUND) return last_vertex;
            }
        }
        colours[v] = BLACK;
        return NOT_FOUND;
    }

    void dfs_find_briges(const Graph::Vertex &v, const Graph::Vertex &ancestor, const Graph &graph,
                         std::vector<vertexColours> &colours, std::set<Edge> &bridges, size_t &timer,
                         std::vector<int> &tin, std::vector<int> &tup) {
        colours[v] = GRAY;
        tup[v] = tin[v] = timer++;
        for (const Graph::Vertex &to : graph.get_all_neighbours(v)) {
            if (to == ancestor) {
                continue;
            } else if (colours[to] == GRAY) {
                tup[v] = std::min(tup[v], tin[to]);
            } else if (colours[to] == WHITE) {
                dfs_find_briges(to, v, graph, colours, bridges, timer, tin, tup);
                tup[v] = std::min(tup[to], tup[v]);
                auto nbrs = graph.get_all_neighbours(v);
                if (tup[to] > tin[v] &&
                    std::count(nbrs.begin(), nbrs.end(), to) == 1) {
                    bridges.insert({v, to});
                }
            }
        }
    }

    void dfs_find_cutpoints(const Graph::Vertex &v, const Graph::Vertex &ancestor, const Graph &graph,
                            std::vector<vertexColours> &colours, std::set<Graph::Vertex> &cutpoints, size_t &timer,
                            std::vector<int> &tin, std::vector<int> &tup) {
        colours[v] = GRAY;
        tup[v] = tin[v] = timer++;
        int children = 0;
        for (const Graph::Vertex &to : graph.get_all_neighbours(v)) {
            if (to == ancestor) {
                continue;
            } else if (colours[to] == GRAY) {
                tup[v] = std::min(tup[v], tin[to]);
            } else if (colours[to] == WHITE) {
                children++;
                dfs_find_cutpoints(to, v, graph, colours, cutpoints, timer, tin, tup);
                tup[v] = std::min(tup[to], tup[v]);
                if (tup[to] >= tin[v] && ancestor != graph.get_vertex_count()) { // is not root
                    cutpoints.insert(v);
                }
            }
        }
        if (ancestor == graph.get_vertex_count() && children > 1) { // is root
            cutpoints.insert(v);
        }
    }


    std::vector<Component> find_components(const Graph &graph) {
        std::vector<Component> connected_components;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                connected_components.emplace_back(std::vector<Graph::Vertex>());
                dfs_find_components(i, graph, colours, connected_components.back());
            }
        }

        return connected_components;
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
        std::vector<Graph::Vertex> result_path;
        auto cur = finish;
        for (int i = 0; i < distance[finish] + 1; i++) {
            result_path.push_back(cur);
            cur = ancestor[cur];
        }
        std::reverse(result_path.begin(), result_path.end());
        return result_path;
    }

    std::vector<Graph::Vertex> graph_top_sort(const Graph &graph) {
        std::vector<Graph::Vertex> top_sort_result;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                dfs_graph_top_sort(i, graph, colours, top_sort_result);
            }
        }

        return top_sort_result;
    }

    std::vector<Graph::Vertex> find_cicle(const Graph &graph) {
        const Graph::Vertex NOT_FOUND = graph.get_vertex_count();
        const Graph::Vertex NOT_SET = graph.get_vertex_count();
        Graph::Vertex last_vertex = NOT_FOUND;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        std::vector<Graph::Vertex> ancestors(graph.get_vertex_count(), NOT_SET);

        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                last_vertex = dfs_find_cicle(i, graph, colours, ancestors);
            }
            if (last_vertex != NOT_FOUND) {
                break;
            }
        }

        std::vector<Graph::Vertex> cicle;
        if (last_vertex == NOT_FOUND) {
            return cicle;
        }
        cicle.push_back(last_vertex);
        auto cur = ancestors[last_vertex];
        while (cur != last_vertex) {
            cicle.push_back(cur);
            cur = ancestors[cur];
        }
        std::reverse(cicle.begin(), cicle.end());
        return cicle;
    }

    std::set<Edge> find_bridges(const Graph &graph) {
        const int NOT_SET = -1;
        const Graph::Vertex NOT_AC = graph.get_vertex_count();
        std::set<Edge> bridges;
        size_t timer = 0;
        std::vector<int> tin(graph.get_vertex_count(), NOT_SET), tup(graph.get_vertex_count(), NOT_SET);
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                dfs_find_briges(i, NOT_AC, graph, colours, bridges, timer, tin, tup);
            }
        }
        return bridges;
    }

    std::set<Graph::Vertex> find_cutpoints(const Graph &graph) {
        const int NOT_SET = -1;
        const Graph::Vertex NOT_AC = graph.get_vertex_count();
        std::set<Graph::Vertex> cutpoints;
        size_t timer = 0;
        std::vector<int> tin(graph.get_vertex_count(), NOT_SET), tup(graph.get_vertex_count(), NOT_SET);
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                dfs_find_cutpoints(i, NOT_AC, graph, colours, cutpoints, timer, tin, tup);
            }
        }
        return cutpoints;
    }

    bool is_bipartite(const Graph &graph) {
        std::vector<int> part(graph.get_vertex_count(), 0);
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);

        for (int i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                part[i] = 1;
                if (!dfs_bipartite(i, graph, colours, part)) {
                    return false;
                }
            }
        }
        return true;
    }
};

template<class T, class Hash = std::hash<T>>
class DisjointSet {
public:
    class DSUExeption {
    };

    DisjointSet(std::vector<T> _v) {
        for (const T &elem : _v) {
            size[elem] = 1;
            prev[elem] = elem;
        }
    }

    DisjointSet(T start, T finish) { // for iterable type, not a reference since copy is needed for cycle
        for (; start <= finish; start++) {
            size[start] = 1;
            prev[start] = start;
        }
    }

    void union_sets(const T &first, const T &second) {
        T first_root = get_set_id(first);
        T second_root = get_set_id(second);
        if (first_root == second_root) {
            throw DSUExeption();
        }
        if (size[second_root] > size[first_root]) {
            // in case std::swap is undefined for T
            T tmp = first_root;
            first_root = second_root;
            second_root = tmp;
        }
        prev[second_root] = first_root;
        size[first_root] += size[second_root];

    }

    size_t get_set_size(const T &elem) {
        return size[get_set_id(elem)];
    }

    bool in_same_set(const T &first, const T &second) {
        return get_set_id(first) == get_set_id(second);
    }

private:
    std::unordered_map<T, size_t, Hash> size;
    std::unordered_map<T, T, Hash> prev;

    T get_set_id(const T &elem) {
        if (elem == prev[elem]) {
            return elem;
        }
        return prev[elem] = get_set_id(prev[elem]);
    }
};


int main() {
    size_t vertex_count, edges_count;
    std::cin >> vertex_count >> edges_count;

    AdjListGraph graph(vertex_count, false);

    for (int i = 0; i < edges_count; i++) {
        Graph::Vertex start, end;
        std::cin >> start >> end;
        graph.add_edge(start - 1, end - 1);
    }

    if (GraphAlgorithms::is_bipartite(graph)) {
        std::cout << "YES\n";
    } else {
        std::cout << "NO\n";
    }
}
