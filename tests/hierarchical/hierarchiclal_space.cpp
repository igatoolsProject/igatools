
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
 * 3. Design of HSpace class and accessor (starting from this file)
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
 * A box is passed and returned as the tensor indices of two opposite corners
 * I will move a list of elements (elems_to_move) that has to be initialized
 *  twice for every direction (moving left and right), and is overwritten at
 *  every step.
 * The difficulties in igatools may be how to move inside the list of
 *  elems_to_move, to store the new element in the right position (in the
 *  pseudo-code the indices [i,j,k]); and how to compute elem_new, moving
 *  with the tensor index ti inside the grid.
 *
 *   MaximalBox(indices)
 *   {
 *    for idim = 1:ndim
 *      for dir = {-1, 1} (left-right, up-down, front-rear, back-forth...)
 *        initialize elems_to_move  (from indices, idim and dir)
 *        intialize ti (from idim and idir, {1,0,0},{-1,0,0},{0,1,0}...)
 *        add_line = true;
 *        while (add_line)
 *          for (elem in elems_to_move)
 *            elem_new = elem++; (elem->move, using ti)
 *            if (!is_influence(elem_new) .or. is_null(elem_new))
 *              add_line = false;
 *              break;
 *            else
 *              elems_to_move[i,j,k] = elem_new;  (overwrite, if possible)
 *            endif
 *          endfor
 *          if (add_line)
 *            indices[idim,max-min] += dir; (max or min, depending on dir)
 *          endif
 *        endwhile
 *      endfor
 *    endfor
 *   }
 *
 * 5. Computation of the active and inactive functions
 *    The algorithm works for any given mesh (either refining or coarsening)
 *
 *    Prerequisite: All maximal boxes from 0 to N should have been computed
 *    N = finest level after the adaptive modification
 *    Set Maxboxes(N+1) = empty;
 *
 *    for l=N:-1:0
 *      - Initialize all the functions to inactive
 *      - Maxboxes_tmp[l+1] = transform index of Maxboxes[l+1] to level l
 *      - Activate functions in V[l] using Maxboxes[l]
 *      - Deactivate functions in V[l] using Maxboxes_tmp[l+1]
 *    end
 *
 */
template<typename Space>
class HierarchicalSpace
{
public:

    using ElementAccessor  = HierarchicalSpaceElementAccessor<self_t>;
    using ElementIterator  = CartesianGridIterator<ElementAccessor>;

    using SpaceElementIterator =  typename Space::ElementIterator;
public:


    /**
     * - activate new elements on the fine mesh / deactivate on the coarse / mark influence (only on the active elements)
     * - activate deactivate functions
     *
     *
     *
     */
    void refine()
    {
        //1. activate, deactivate and mark influence on elements;
        //2. compute maximal boxes
        //3. activate/deactivate function
    }

private:
    SpaceElementIterator parent(const SpaceElementIterator &elem);

private:
    Tree<dim> refinement_tree_;

    vector<Space> spaces_;
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
    }
    while (elem != NULL);

    return res;
}
