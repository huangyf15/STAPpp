*Heading
** some comment

*Preprint, echo=NO, model=NO, history=NO, contact=NO

*Part, name=$NAME
*End Part

*Assembly, name=$NAME

    *Instance, name=$NAME, part=$name
        *Node
            {nodeIndex,  x, y, z}
        *Element, type=$TYPE
            {elementIndex, {indexOfNodesOfElement}}

        *Nset, nset=$NAME, internal, generate
            from, to, step
        *Elset, elset=$NAME, internal, generate
            from, to, step
        *Solid Section, elset=$ELSET_NAME, material=$MATERIAL_NAME
        ,
    *End Instance

    *Nset, nset=$NAME, internal, instance=$NAME
        {nodeIndex}
    *Elset, elset=$NAME, internal, instance=$NAME, generate
        from, to, step
    *Surface, type=$TYPE, name=$NAME, internal
        $SET_NAME, $DIRECTION

*End Assembly

*Material, name=$NAME
*Elastic
    {E, poison rate}

**------------------------------------------------------

*Step, name=$NAME, nlgeom=$YN
    *Static
        $INIT, $TIME_PERIOD, $MIN, $MAX

    *Boundary
        $SET_NAME, $DOF, $DOF

    *DsLoad
        $SET_NAME, $CODE, $MAG
        ** _PickedSurf5, P, 1.

    *Restart, write, frequency=0

    *Output, field, variable=$NAME
    *Output, history, variable=$NAME
*End Step

