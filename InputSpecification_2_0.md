# Input version 2.0 specification #

The igatools input files are defined using the Extensible Markup Language 1.0 ([XML](http://www.w3.org/TR/REC-xml/,)).

Below, an example of igatools input (version 2.0) of an `IgMapping` is shown. The `IgMapping` XML tag is splitted in two main parts:
  * Reference space: the `NURBSSpace` tag
  * The control points: `ControlPoints` tag


```
<?xml version="1.0" encoding="utf-8"?>
<Igatools FormatVersion="2.0">
  <IgMapping Dim="2" Codim="0" RefSpaceType="NURBSSpace">
    <NURBSSpace Dim="2" Range="2" Rank="1">
      <CartesianGrid Dim="2">
        <Knots Direction="0" Size="3">
          0.000000 0.500000 1.000000
        </Knots>
        <Knots Direction="1" Size="3">
          0.000000 0.500000 1.000000
        </Knots>
      </CartesianGrid>
      <NURBSSpaceScalarComponents Size="2">
        <ComponentsMap Size="2">
          0 1
        </ComponentsMap>
        <NURBSSpaceScalarComponent Id="0">
          <DofsTensorSize Dim="2">
            4 4
          </DofsTensorSize>
          <Degrees Dim="2">
            2 2
          </Degrees>
          <InteriorMultiplicities Dim="2">
            <InteriorMultiplicity Direction="0" Size="1">
              1
            </InteriorMultiplicity>
            <InteriorMultiplicity Direction="1" Size="1">
              1
            </InteriorMultiplicity>
          </InteriorMultiplicities>
          <Weights Size="16">
            1.000000 1.000000 1.000000 1.000000
            0.853553 0.853553 0.853553 0.853553
            0.853553 0.853553 0.853553 0.853553
            1.000000 1.000000 1.000000 1.000000
          </Weights>
        </NURBSSpaceScalarComponent>
        <NURBSSpaceScalarComponent Id="1">
          <DofsTensorSize Dim="2">
            4 4
          </DofsTensorSize>
          <Degrees Dim="2">
            2 2
          </Degrees>
          <InteriorMultiplicities Dim="2">
            <InteriorMultiplicity Direction="0" Size="1">
              1
            </InteriorMultiplicity>
            <InteriorMultiplicity Direction="1" Size="1">
              1
            </InteriorMultiplicity>
          </InteriorMultiplicities>
          <Weights Size="16">
            1.000000 1.000000 1.000000 1.000000
            0.853553 0.853553 0.853553 0.853553
            0.853553 0.853553 0.853553 0.853553
            1.000000 1.000000 1.000000 1.000000
          </Weights>
        </NURBSSpaceScalarComponent>
      </NURBSSpaceScalarComponents>
    </NURBSSpace>
    <ControlPoints Dim="1" Size="32">
      1.000000 1.750000 3.250000 4.000000
      1.000000 1.750000 3.250000 4.000000
      0.414214 0.724874 1.346194 1.656854
      0.000000 0.000000 0.000000 0.000000
      0.000000 0.000000 0.000000 0.000000
      0.414214 0.724874 1.346194 1.656854
      1.000000 1.750000 3.250000 4.000000
      1.000000 1.750000 3.250000 4.000000
    </ControlPoints>
  </IgMapping>
</Igatools>
```

The example above is mostly self-explanatory but some details are
worthy of comment:

  * The values of `ControlPoints` are written component by component, i.e. first all the values for the first component are written, then the second component, and so on.
  * It is possible to create reference spaces with more than one component without defining all of them explicitly, specifing them in the `ComponentsMap` tag. Let's see in an example:
```
<BSplineSpace Dim="3" Range="3" Rank="1">
  <CartesianGrid Dim="3">
    <Knots Direction="0" Size="3">
      0.000000 0.500000 1.000000
    </Knots>
    <Knots Direction="1" Size="3">
      0.000000 0.500000 1.000000
    </Knots>
    <Knots Direction="2" Size="4">
      0.000000 0.333333 0.6666660 1.000000
    </Knots>
  </CartesianGrid>
  <BSplineSpaceScalarComponents Size="1">
    <ComponentsMap Size="3">
      0 0 0
    </ComponentsMap>
    <BSplineSpaceScalarComponent Id="0">
      <DofsTensorSize Dim="3">
        4 4
      </DofsTensorSize>
      <Degrees Dim="3">
        2 2 3
      </Degrees>
      <InteriorMultiplicities Dim="3">
        <InteriorMultiplicity Direction="0" Size="1">
          1
        </InteriorMultiplicity>
        <InteriorMultiplicity Direction="1" Size="1">
          1
        </InteriorMultiplicity>
        <InteriorMultiplicity Direction="2" Size="2">
          1 2
        </InteriorMultiplicity>
      </InteriorMultiplicities>
    </BSplineSpaceScalarComponent>
  </BSplineSpaceScalarComponents>
</BSplineSpace>
```

In the example above, a `BSplineSpace` with `Range="3"` is built
by defining a single `BSplineSpaceScalarComponent`.
In `ComponentsMap` is specified that the defined component `0` will be used for the three components of the space.