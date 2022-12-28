module octree_mod

  implicit none

  private

  public point_type
  public octree_init
  public octree_final
  public octree_build
  public octree_update
  public octree_search

  type config_type
    integer max_num_point ! minimum number of points we want to be contained in leaf node.
    integer max_depth     ! Maximum level of branch and leaf nodes
    real(8) bbox(2, 3)    ! The 2 points in a box one being the anchor point bottom left then 2nd being its adjacent corner
  end type config_type

  ! Points should be indexed by their id.
  type point_type
    integer id  ! (Assuming) this is to determine which order the numbers are in
    real(8) x(3)  !point is in form (x, y, z)
  end type point_type

  ! There are two kinds of nodes:
  !   1. Branch node with children;
  !   2. Leaf node without child but containing points.
  type node_type
    integer depth ! depth tells how far we are into the tree 
    real(8) bbox(2, 3) ! Same as before 2 coordinates of box in form (x,y,z)
    integer num_point ! used to determine the number of points (inside a box? Not quite sure)
    integer, allocatable :: point_ids(:) ! id for the point we are working with
    type(node_type), pointer :: parent  ! Pointer for parent node
    type(node_type), pointer :: children(:) ! points to the 8 children after being created
  end type node_type  ! Change to add morton code possibly to find the cell in one big array.

  type tree_type
    type(point_type), pointer :: points(:)  ! A pointer that points to the points within that cube in the octree
    type(node_type), pointer :: root_node   ! A pointer that points to the root that we are on
  end type tree_type

  type(config_type) config ! Type of configuration to set max num of points and depth  
  type(tree_type) tree     ! Type of configuration to set this tree point and node types for new Children or branches

contains

  subroutine octree_init(max_num_point, max_depth, bbox)  ! Subroutine so that we create the first octree box?

    integer, intent(in), optional :: max_num_point 
    integer, intent(in), optional :: max_depth
    real(8), intent(in), optional :: bbox(2, 3)

    config%max_num_point = merge(max_num_point, 3, present(max_num_point)) ! sets the max number of points if given in the argument otherwise it is set to 3
    config%max_depth = merge(max_depth, 10, present(max_depth)) ! sets the max depth of our octree if given otherwise it is set to 10
    config%bbox = merge(bbox, reshape([0.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0], [2, 3]), present(bbox)) ! creates the first octree box if given the 2 corners with bbox otherwise does a generic box of 0x1 in all 3 dimensions

    if (.not. associated(tree%root_node)) allocate(tree%root_node) ! determines if the tree%root_node has a pointer if false it allocates it
    call reset_node(tree%root_node) ! See Subroutine for details
    tree%root_node%depth = 1 ! Set the root node depth to 1 since we are creating the first tree
    tree%root_node%bbox = config%bbox ! Sets the first node with the boundary box value to determine size of box.

  end subroutine octree_init

  subroutine octree_final()

    call clean_node(tree%root_node) ! See clean_node for details
    deallocate(tree%root_node) ! deallocates tree root node such that we can tell our refinement is over.

  end subroutine octree_final

  recursive subroutine octree_build(points, node_) ! Recursive subroutine that builds the octree for each refinement

    type(point_type), intent(in), target :: points(:) ! This argument is needed when calling the octree_build these are the points in the octree mesh
    type(node_type), intent(inout), target, optional :: node_ ! optional argument to determine if we are in a refine octree cell.

    type(node_type), pointer :: node ! create an empty node pointer
    integer i, j !dummy integers for the recursive subroutine
    integer num_contained_point ! value for number of points contained in our bbox may be bigger than max points
    type(point_type), allocatable :: contained_points(:) !array of the points that was contained in out octree.

    if (present(node_)) then  ! Sets the empty node pointer to point to the argument one in subroutine
      node => node_
    else
      tree%points => points !if the optional node_ is not in the argument, we point the tree points to our argument points  
      node => tree%root_node  ! then the node pointer points to our tree root nodes. 
    end if

    ! Leaf node is approached.
    if (node%depth >= config%max_depth .or. size(points) <= config%max_num_point) then  ! Checks depth node to max depth and size of points, if one is true goes to next if statement
      if (size(points) > size(node%point_ids)) then ! if num of points is more than the number of ids we de allocate the ids and make new array with correct size
        deallocate(node%point_ids)
        allocate(node%point_ids(size(points)))
      end if
      j = 1 ! first if is false for both conditions goes on to this do statement
      do i = 1, size(points) 
        if (points(i)%x(1) < node%bbox(1, 1) .or. points(i)%x(1) > node%bbox(2, 1) .or. & ! checks if the point is inside our bbox in all 3 dimensions
            points(i)%x(2) < node%bbox(1, 2) .or. points(i)%x(2) > node%bbox(2, 2) .or. & ! if the point is outside of box from one side then we cycle the i integer
            points(i)%x(3) < node%bbox(1, 3) .or. points(i)%x(3) > node%bbox(2, 3)) cycle
        node%point_ids(j) = points(i)%id  ! once we find a point inside our box then we set the node id with this points id so we can see what points we retain.
        j = j + 1
      end do
      node%num_point = j - 1 ! sets the total number of points contained in our node.
      return
    end if

    ! Copy contained points into a new array.
    num_contained_point = 0 ! sets a dummy variable to 0 for number of points in our nodes
    do i = 1, size(points)
      if (points(i)%x(1) < node%bbox(1, 1) .or. points(i)%x(1) > node%bbox(2, 1) .or. & ! checks again for the number of points within the cell 
          points(i)%x(2) < node%bbox(1, 2) .or. points(i)%x(2) > node%bbox(2, 2) .or. &
          points(i)%x(3) < node%bbox(1, 3) .or. points(i)%x(3) > node%bbox(2, 3)) cycle
      num_contained_point = num_contained_point + 1 !adds 1 to every point that is contained in our cell.
    end do
    allocate(contained_points(num_contained_point))
    j = 1
    do i = 1, size(points)
      if (points(i)%x(1) < node%bbox(1, 1) .or. points(i)%x(1) > node%bbox(2, 1) .or. & !3rd check where in this loop it assigns the contained points to an array with their ids
          points(i)%x(2) < node%bbox(1, 2) .or. points(i)%x(2) > node%bbox(2, 2) .or. & ! also stores the points coordinates so we can see which cell they lie in.
          points(i)%x(3) < node%bbox(1, 3) .or. points(i)%x(3) > node%bbox(2, 3)) cycle
      contained_points(j)%id = points(i)%id
      contained_points(j)%x = points(i)%x
      j = j + 1
    end do

    if (num_contained_point == 0) return ! on more refinements once we have no contained points in our cells it ends the loop

    ! Subdivide node and run into the child nodes.
    call subdivide_node(node) ! See Subroutine for details, essentially refines to 8 children cells
    do i = 1, 8
      call octree_build(contained_points, node%children(i)) ! builds our octree again now with more refined cells 
    end do

    ! if (node%depth == 1) then
    !   call print_tree(tree%root_node)
    ! end if

  end subroutine octree_build

  subroutine octree_update(node_)

    type(node_type), intent(inout), target, optional :: node_

    type(node_type), pointer :: node

    if (present(node_)) then
      node => node_
    else
      node => tree%root_node
    end if

  end subroutine octree_update

  recursive subroutine octree_search(x, distance, num_ngb_point, ngb_ids, node_)

    real(8), intent(in) :: x(3)
    real(8), intent(in) :: distance
    integer, intent(inout) :: num_ngb_point
    integer, intent(inout) :: ngb_ids(:)
    type(node_type), intent(in), target, optional :: node_

    type(node_type), pointer :: node
    real(8) d2, dx(3)
    integer i

    if (present(node_)) then
      node => node_
    else
      node => tree%root_node
    end if

    if (associated(node%children)) then
      ! We are at branch node.
      do i = 1, 8
        if ((x(1) + distance) > node%children(i)%bbox(1,1) .and. &
            (x(1) - distance) < node%children(i)%bbox(2,1) .and. &
            (x(2) + distance) > node%children(i)%bbox(1,2) .and. &
            (x(2) - distance) < node%children(i)%bbox(2,2) .and. &
            (x(3) + distance) > node%children(i)%bbox(1,3) .and. &
            (x(3) - distance) < node%children(i)%bbox(2,3)) then
          call octree_search(x, distance, num_ngb_point, ngb_ids, node%children(i))
        end if
      end do
    else
      if (node%num_point == 0) return
      ! We are at leaf node.
      d2 = distance * distance
      do i = 1, node%num_point
        dx(:) = x(:) - tree%points(node%point_ids(i))%x(:)
        if (dot_product(dx, dx) < d2) then
          num_ngb_point = num_ngb_point + 1
          if (num_ngb_point <= size(ngb_ids)) then
            ngb_ids(num_ngb_point) = node%point_ids(i)
          else
            write(6, "('[Error]: octree: The ngb_ids array size is not enough!')")
            stop 1
          end if
        end if
      end do
    end if

  end subroutine octree_search

  subroutine reset_node(node)

    type(node_type), intent(inout) :: node

    node%num_point = 0  ! sets number of points in node to 0
    if (.not. allocated(node%point_ids)) allocate(node%point_ids(config%max_num_point)) ! Checks if the node point ids is allocated, if it isnt it allocates it with the minimum number of points in our cell we wish to have
    nullify(node%parent)  ! takes away pointer of node to the parents
    if (size(node%children) == 8) deallocate(node%children) ! Checks if we already have children cells, if true we deallocate the node pointer for the children
    nullify(node%children)  ! takes away pointer of node for children

  end subroutine reset_node

  subroutine subdivide_node(node)

    type(node_type), intent(inout), target :: node ! this argument is not optional and subroutine is called once we have refined our octree atleast once.

    integer i, j, k, l  ! creates 4 interger so we can run through all children cells and to check each dimension of the cell that we are on.
    real(8) bbox(2, 3)

    allocate(node%children(8)) ! fills the children pointer with empty array.
    l = 1 ! First Child cell, and continues until all cells are approached
    do k = 1, 2 ! z dimension
      do j = 1, 2 ! y dimension
        do i = 1, 2 ! x dimension
          call reset_node(node%children(l)) ! runs subroutine for reset_node directly above (essentially resets pointers and deallocates children nodes
          node%children(l)%depth = node%depth + 1 ! gives the depth for children l cell to be next level of depth
          node%children(l)%parent => node ! re assigns parent pointer with new children cell towards the node argument we input.  
          node%children(l)%bbox(1, 1) = node%bbox(1, 1) + (i - 1) * (node%bbox(2, 1) - node%bbox(1, 1)) * 0.5d0 ! The following lines divide each dimension in half so that we add 
          node%children(l)%bbox(2, 1) = node%bbox(2, 1) - (2 - i) * (node%bbox(2, 1) - node%bbox(1, 1)) * 0.5d0 ! a total of 8 cells each of them having their own respective pointers
          node%children(l)%bbox(1, 2) = node%bbox(1, 2) + (j - 1) * (node%bbox(2, 2) - node%bbox(1, 2)) * 0.5d0 ! along with become new "parent" cells if further refined
          node%children(l)%bbox(2, 2) = node%bbox(2, 2) - (2 - j) * (node%bbox(2, 2) - node%bbox(1, 2)) * 0.5d0 ! it gives new coordinates for each anchor point and diagonal point
          node%children(l)%bbox(1, 3) = node%bbox(1, 3) + (k - 1) * (node%bbox(2, 3) - node%bbox(1, 3)) * 0.5d0 ! in order to create our new refined cells.
          node%children(l)%bbox(2, 3) = node%bbox(2, 3) - (2 - k) * (node%bbox(2, 3) - node%bbox(1, 3)) * 0.5d0
          node%children(l)%parent => node ! same? not sure why we need it twice.
          l = l + 1 ! increases child cell number until we hit 8
        end do
      end do
    end do

  end subroutine subdivide_node

  recursive subroutine clean_node(node)

    type(node_type), intent(inout) :: node

    integer i

    if (associated(node%children)) then
      do i = 1, 8
        call clean_node(node%children(i))
        deallocate(node%children(i)%point_ids)
      end do
      deallocate(node%children)
    end if

  end subroutine clean_node

  subroutine print_node(node)

    type(node_type), intent(in) :: node

    write(6, "('Bounding box: ', 6F8.2)") node%bbox
    write(6, "('Depth: ', I3)") node%depth
    write(6, "('Point number: ', I3)") node%num_point
    write(6, "('Leaf?: ', L1)") .not. associated(node%children)

  end subroutine print_node

  recursive subroutine print_tree(node)

    type(node_type), intent(in) :: node

    integer i

    if (associated(node%children)) then
      write(6, "('----------------------------------------------------------------')")
      write(6, "('Branch node: ')")
      write(6, "('  Bounding box: ', 6F8.2)") node%bbox
      write(6, "('  Depth: ', I3)") node%depth
      do i = 1, 8
        call print_tree(node%children(i))
      end do
    else
      if (node%num_point == 0) return
      write(6, "('Leaf node: ')")
      write(6, "('  Bounding box: ', 6F8.2)") node%bbox
      write(6, "('  Depth: ', I3)") node%depth
      write(6, "('  Points:')", advance='no')
      write(6, *) (node%point_ids(i), i = 1, node%num_point)
    end if

  end subroutine print_tree

end module octree_mod
