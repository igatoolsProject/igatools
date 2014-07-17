
/**
 *
 * 1. Make CartesianGrid aware of:
 *   a)
 *     - active elements
 *     - influence elements
 *     - marked elements (different types)
 *
 *   b)
 *     - get_num_elem() = n. active elements
 *     - get_num_all_elements() = n. active + influence
 *
 *   c) CartesianGridElementAccessor:
 *     - Modify operator++ to increment by active elements
 *     - check index access
 *
 * 2. Make Spaces aware of active functions
 *   a)
 *     - vector of active functions
 *     - get_num_functions() = n. active functions
 *  
 *   b) Accessor, all local number of functions should be revisit
 *     - get_local_to_global
 *     - get_values
 *
 * 3. Design of HSpace class and accessor
 *
 * 4. Design of the MaximalBox algorithm
 *  a) Given a CartesianGrid with the influence elements information, 
 *     find a maximal box covering: 
 *
 *     list_yet_to_cover = influence elements on CartesianGrid
 *     Repeat:
 *       - \Omega_e \in list
 *       - Box = MaximalBox(\Omega_e)
 *       - list = list - Box
 *     until(list == empty)
 *
 *     MaximalBox(\Omega_e)
 *     {
 *       Move L and R(\Omega_e) -> elem_list_0
 *       Move U and D elem_list_0 -> elem_list_1
 *       .....
 *       .....
 *       return a Box (i.e. the tensor indices of the opposite 
 *       elements at the corners defining the Box)
 *     }
 *
 *
 *
 *
 *
 */
template<typename Space>
class HierarchicalSpace
{
public:

  using ElementAccessor  = HierarchicalSpaceElementAccessor<self_t>;
  using ElementIterator  = GridForwardIterator<ElementAccessor>;

  using SpaceElementIterator =  typename Space::ElementIterator;
public:
  

  /**
   * - activate new elements on the fine mesh / deactivate on the coarse / mark influence (only on the active elements)
   * - activate deactivate functions
   *
   *
   *
   */
  void refine();

private:
  SpaceElementIterator parent(const SpaceElementIterator &elem);

private:
  Tree<dim> refinement_tree_;
  
  std::vector<Space> spaces_;
}


template<typename Space>
class HierarchicalSpaceAccessor : public typename Space::Accessor;
{

public:

  ValueTable<Value> 
  get_basis_values(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;
  
private:
  const std::shared_ptr<ContainerType> h_space_;
  int level_;
  SpaceIterator elem_;
}


auto
HierarchicalSpaceAccessor<Space>::
get_basis_values(const TopologyId<dim> &topology_id) const -> ValueTable<Value> 
{
  ValueTable<Value> res;
  auto elem = elem_;
  do {
    res += elem->get_basis_values();
    elem = h_space_->parent(elem);
  } while (elem != NULL);

  return res;
}
