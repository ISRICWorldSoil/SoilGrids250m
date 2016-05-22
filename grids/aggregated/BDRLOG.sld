<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
  <NamedLayer>
    <Name>SoilGrids250m</Name>
    <UserStyle>
      <Title>BDRLOG_M_250m</Title>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
            <Geometry>
              <PropertyName>Predicted probability of occurence (0-100%) of R horizon</PropertyName>
              <Opacity>1</Opacity>
            </Geometry>
            <ColorMap type="intervals">
              <ColorMapEntry color="#FFFFFF" label="NODATA" opacity="0.0" quantity="255.0"/>
              <ColorMapEntry color="#FFFFCC" label="0" opacity="0.7" quantity="0.0"/>
              <ColorMapEntry color="#FFF0A8" label="10" opacity="0.7" quantity="10.0"/>
              <ColorMapEntry color="#FEE186" label="20" opacity="0.7" quantity="20.0"/>
              <ColorMapEntry color="#FEC965" label="30" opacity="0.7" quantity="30.0"/>
              <ColorMapEntry color="#FDAA48" label="40" opacity="0.7" quantity="40.0"/>
              <ColorMapEntry color="#FD8D3C" label="50" opacity="0.7" quantity="50.0"/>
              <ColorMapEntry color="#FC5A2D" label="60" opacity="0.7" quantity="60.0"/>
              <ColorMapEntry color="#EC2E21" label="70" opacity="0.7" quantity="70.0"/>
              <ColorMapEntry color="#D30F20" label="80" opacity="0.7" quantity="80.0"/>
              <ColorMapEntry color="#B00026" label="90" opacity="0.7" quantity="90.0"/>
              <ColorMapEntry color="#800026" label="100" opacity="0.7" quantity="100.0"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>