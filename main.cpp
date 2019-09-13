#include <iostream>
#include <vector>

class Graph {
public:
    typedef size_t Vertex;
    Graph(size_t _vertex_count, bool _is_directed) {
        is_directed = _is_directed;
        vertex_count = _vertex_count;
    }

    bool getIsDirected() {
        return is_directed;
    }
    size_t getVertexCount() {
        return vertex_count;
    }
    size_t getEdgesCount() {
        return edges_count;
    }

    virtual void addEdge(const Vertex & start, const Vertex & finish) {
        edges_count++;
    }
protected:
    size_t vertex_count, edges_count;
    bool is_directed;

    virtual std::vector<Vertex> getAllNeighbours(Vertex v) const = 0;

};
class AdjListGraph: public Graph {
public:
    AdjListGraph(size_t _vertex_count, bool _is_directed):Graph(_vertex_count, _is_directed) {
        Vertexes.resize(vertex_count);
    }
    virtual void addEdge(const Vertex & start, const Vertex & finish) override {
        Vertexes[start].push_back(finish);
        if(!is_directed) {
            Vertexes[finish].push_back(start);
        }
    }
    size_t GetVertexDeg(Vertex v) {
        return Vertexes[v].size();
    }
private:
    virtual std::vector<Vertex> getAllNeighbours(Vertex v) const override{
        return std::vector<Vertex>();
    }
    std::vector<std::vector<Vertex>> Vertexes;
};


int main() {
    int n, m;
    std::cin >> n >> m;
    AdjListGraph v(n, false);
    for (int i = 0; i < m; i++) {
        Graph::Vertex s, f;
        std::cin >> s >> f;
        v.addEdge(s - 1, f - 1);
    }

    for (int i = 0; i < n; i++) {
        std::cout << v.GetVertexDeg(i) << std::endl;
    }
    return 0;
}

