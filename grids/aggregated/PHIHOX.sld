<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
  <NamedLayer>
    <Name>SoilGrids250m</Name>
    <UserStyle>
      <Title>PHIHOX</Title>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
           <Geometry>
              <PropertyName>Soil pH x 10 in H2O</PropertyName>
              <Opacity>1</Opacity>
            </Geometry>
            <ColorMap type="intervals">
              <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
              <ColorMapEntry color="#FF0000" quantity="42" label="20 - 42" opacity="0.7"/>
              <ColorMapEntry color="#FF1C00" quantity="45" label="42 - 45" opacity="0.7"/>
              <ColorMapEntry color="#FF3900" quantity="46" label="45 - 46" opacity="0.7"/>
              <ColorMapEntry color="#FF5500" quantity="48" label="46 - 48" opacity="0.7"/>
              <ColorMapEntry color="#FF7100" quantity="49" label="48 - 49" opacity="0.7"/>
              <ColorMapEntry color="#FF8E00" quantity="50" label="49 - 50" opacity="0.7"/>
              <ColorMapEntry color="#FFAA00" quantity="51" label="50 - 51" opacity="0.7"/>
              <ColorMapEntry color="#FFC600" quantity="52" label="51 - 52" opacity="0.7"/>
              <ColorMapEntry color="#FFE200" quantity="53" label="52 - 53" opacity="0.7"/>
              <ColorMapEntry color="#FFFF00" quantity="54" label="53 - 54" opacity="0.7"/>
              <ColorMapEntry color="#E3FF00" quantity="55" label="54 - 55" opacity="0.7"/>
              <ColorMapEntry color="#C7FF00" quantity="56" label="55 - 56" opacity="0.7"/>
              <ColorMapEntry color="#AAFF00" quantity="57" label="56 - 57" opacity="0.7"/>
              <ColorMapEntry color="#8EFF00" quantity="58" label="57 - 58" opacity="0.7"/>
              <ColorMapEntry color="#72FF00" quantity="59" label="58 - 59" opacity="0.7"/>
              <ColorMapEntry color="#55FF00" quantity="60" label="59 - 60" opacity="0.7"/>
              <ColorMapEntry color="#39FF00" quantity="61" label="60 - 61" opacity="0.7"/>
              <ColorMapEntry color="#1DFF00" quantity="62" label="61 - 62" opacity="0.7"/>
              <ColorMapEntry color="#01FF00" quantity="63" label="62 - 63" opacity="0.7"/>
              <ColorMapEntry color="#00FF1C" quantity="64" label="63 - 64" opacity="0.7"/>
              <ColorMapEntry color="#00FF38" quantity="65" label="64 - 65" opacity="0.7"/>
              <ColorMapEntry color="#00FF54" quantity="66" label="65 - 66" opacity="0.7"/>
              <ColorMapEntry color="#00FF71" quantity="68" label="66 - 68" opacity="0.7"/>
              <ColorMapEntry color="#00FF8D" quantity="69" label="68 - 69" opacity="0.7"/>
              <ColorMapEntry color="#00FFA9" quantity="70" label="69 - 70" opacity="0.7"/>
              <ColorMapEntry color="#00FFC6" quantity="71" label="70 - 71" opacity="0.7"/>
              <ColorMapEntry color="#00FFE2" quantity="73" label="71 - 73" opacity="0.7"/>
              <ColorMapEntry color="#00FFFE" quantity="75" label="73 - 75" opacity="0.7"/>
              <ColorMapEntry color="#00E3FF" quantity="76" label="75 - 76" opacity="0.7"/>
              <ColorMapEntry color="#00C7FF" quantity="77" label="76 - 77" opacity="0.7"/>
              <ColorMapEntry color="#00ABFF" quantity="78" label="77 - 78" opacity="0.7"/>
              <ColorMapEntry color="#008FFF" quantity="79" label="78 - 79" opacity="0.7"/>
              <ColorMapEntry color="#0072FF" quantity="80" label="79 - 80" opacity="0.7"/>
              <ColorMapEntry color="#0056FF" quantity="81" label="80 - 81" opacity="0.7"/>
              <ColorMapEntry color="#003AFF" quantity="82" label="81 - 82" opacity="0.7"/>
              <ColorMapEntry color="#001DFF" quantity="83" label="82 - 83" opacity="0.7"/>
              <ColorMapEntry color="#0001FF" quantity="84" label="83 - 84" opacity="0.7"/>
              <ColorMapEntry color="#1B00FF" quantity="86" label="84 - 86" opacity="0.7"/>
              <ColorMapEntry color="#3800FF" quantity="89" label="86 - 89" opacity="0.7"/>
              <ColorMapEntry color="#5400FF" quantity="110" label="89 - 110" opacity="0.7"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>