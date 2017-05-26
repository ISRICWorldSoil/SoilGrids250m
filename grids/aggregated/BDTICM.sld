<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
  <NamedLayer>
    <Name>SoilGrids250m</Name>
    <UserStyle>
      <Title>BDTICM_M_250m</Title>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
            <Geometry>
              <ogc:PropertyName>Absolute depth to bedrock (in cm)</ogc:PropertyName>
            </Geometry>
            <ColorMap type="intervals">
              <ColorMapEntry color="#FFFFFF" label="NODATA" opacity="0.0" quantity="-99999.0"/>
              <ColorMapEntry color="#FFFFD9" label="0" opacity="1.0" quantity="0.0"/>
              <ColorMapEntry color="#F7FCC8" label="30" opacity="1.0" quantity="30.0"/>
              <ColorMapEntry color="#EFF9B7" label="90" opacity="1.0" quantity="90.0"/>
              <ColorMapEntry color="#E3F4B1" label="150" opacity="1.0" quantity="150.0"/>
              <ColorMapEntry color="#D3EDB3" label="240" opacity="1.0" quantity="240.0"/>
              <ColorMapEntry color="#BFE6B4" label="360" opacity="1.0" quantity="360.0"/>
              <ColorMapEntry color="#A1DAB7" label="580" opacity="1.0" quantity="580.0"/>
              <ColorMapEntry color="#82CEBA" label="880" opacity="1.0" quantity="880.0"/>
              <ColorMapEntry color="#68C4BE" label="1050" opacity="1.0" quantity="1050.0"/>
              <ColorMapEntry color="#4EBAC2" label="1200" opacity="1.0" quantity="1200.0"/>
              <ColorMapEntry color="#39AEC3" label="1350" opacity="1.0" quantity="1350.0"/>
              <ColorMapEntry color="#2A9EC1" label="1650" opacity="1.0" quantity="1650.0"/>
              <ColorMapEntry color="#1D8EBE" label="1900" opacity="1.0" quantity="1900.0"/>
              <ColorMapEntry color="#1F78B4" label="2100" opacity="1.0" quantity="2100.0"/>
              <ColorMapEntry color="#2163AA" label="2400" opacity="1.0" quantity="2400.0"/>
              <ColorMapEntry color="#2250A1" label="3200" opacity="1.0" quantity="3200.0"/>
              <ColorMapEntry color="#243F99" label="5400" opacity="1.0" quantity="5400.0"/>
              <ColorMapEntry color="#20308A" label="9100" opacity="1.0" quantity="9100.0"/>
              <ColorMapEntry color="#142671" label="15000" opacity="1.0" quantity="15000.0"/>
              <ColorMapEntry color="#081D58" label="54121" opacity="1.0" quantity="54121.0"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>
