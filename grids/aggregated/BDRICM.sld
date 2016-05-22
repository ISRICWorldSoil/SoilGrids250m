<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
  <NamedLayer>
    <Name>SoilGrids250m</Name>
    <UserStyle>
      <Title>BDRICM_M_250m</Title>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
            <Geometry>
              <PropertyName>Depth to bedrock (R horizon) up to 200 cm</PropertyName>
              <Opacity>1</Opacity>
            </Geometry>
            <ColorMap type="intervals">
              <ColorMapEntry color="#FFFFFF" label="NODATA" opacity="0.0" quantity="255.0"/>
              <ColorMapEntry color="#7F2704" label="0" opacity="0.7" quantity="0.0"/>
              <ColorMapEntry color="#AD3802" label="30" opacity="0.7" quantity="30.0"/>
              <ColorMapEntry color="#DF5106" label="60" opacity="0.7" quantity="60.0"/>
              <ColorMapEntry color="#F67824" label="90" opacity="0.7" quantity="90.0"/>
              <ColorMapEntry color="#FD9F56" label="120" opacity="0.7" quantity="120.0"/>
              <ColorMapEntry color="#FDC692" label="160" opacity="0.7" quantity="160.0"/>
              <ColorMapEntry color="#FDE2C7" label="190" opacity="0.7" quantity="190.0"/>
              <ColorMapEntry color="#FFF5EB" label="200" opacity="0.7" quantity="200.0"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>