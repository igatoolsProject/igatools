template<typename Space>
class HierarchicalSpace
{
public:

  using ElementAccessor  = HierarchicalSpaceElementAccessor<self_t>;
  using ElementIterator  = GridForwardIterator<ElementAccessor>;

  using SpaceElementIterator =  typename Space::ElementIterator;
public:
  
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
