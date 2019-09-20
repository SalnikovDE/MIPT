#include <iostream>
#include <vector>


/*
В галактике «Milky Way» на планете «Snowflake» есть N городов,
 некоторые из которых соединены дорогами. Император галактики «Milky Way» решил
 провести инвентаризацию дорог на планете «Snowflake». Но, как оказалось,
 он не силен в математике, поэтому он просит вас сосчитать количество дорог.
 Требуется написать программу, помогающую императору сосчитать количество дорог на планете «Snowflake».
*/

class Graph {
public:
    typedef size_t Vertex;
    Graph(size_t _vertex_count, bool _is_directed) {
        is_directed = _is_directed;
        vertex_count = _vertex_count;
        edges_count = 0;
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

    virtual std::vector<Vertex> getAllNeighbours(const Vertex &v) const = 0;

};
class AdjListGraph: public Graph {
public:
    AdjListGraph(size_t _vertex_count, bool _is_directed):Graph(_vertex_count, _is_directed) {
        Vertexes.resize(vertex_count);
    }
    virtual void addEdge(const Vertex & start, const Vertex & finish) override {
        Graph::addEdge(start, finish);
        Vertexes[start].push_back(finish);
        if(!is_directed) {
            Vertexes[finish].push_back(start);
        }
    }

    size_t GetVertexDeg(Vertex v) {
        return Vertexes[v].size();
    }

    void ReadMatrix (const std::vector<std::vector<int>> &matrix) {
        for (Vertex i = 0; i < matrix.size(); i++){
            for (Vertex j = 0; j < matrix[i].size(); j++) {
                if (matrix[i][j] != 0) {
                    if(j >= i || is_directed) {
                        addEdge(i, j);
                    }
                }
            }
        }
    };

private:
    virtual std::vector<Vertex> getAllNeighbours(const Vertex &v) const override{
        return std::vector<Vertex>();
    }
    std::vector<std::vector<Vertex>> Vertexes;
};


int main() {
    int n;
    std::cin >> n;
    AdjListGraph g(n, false);
    std::vector<std::vector<int>> v;
    for (int i = 0; i < n; i++) {
        v.emplace_back(std::vector<int>());
        for (int j = 0; j < n; j++) {
            int a;
            std::cin >> a;
            v[i].push_back(a);
        }
    }
    g.ReadMatrix(v);
    std::cout << g.getEdgesCount();
    return 0;
}
