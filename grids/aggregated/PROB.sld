<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
    <UserLayer>
        <LayerFeatureConstraints>
            <FeatureTypeConstraint/>
        </LayerFeatureConstraints>
        <UserStyle>
            <Name>Probability</Name>
            <Title/>
            <FeatureTypeStyle>
                <Name/>
                <Rule>
                    <RasterSymbolizer>
                        <Geometry>
                            <ogc:PropertyName>Predicted probability in percent</ogc:PropertyName>
                        </Geometry>
                        <ColorMap>
                            <ColorMapEntry color="#0000ff" label="0.0" opacity="1.0" quantity="0"/>
                            <ColorMapEntry color="#0028d7" label="2.1" opacity="1.0" quantity="2.10526"/>
                            <ColorMapEntry color="#0050af" label="4.2" opacity="1.0" quantity="4.21053"/>
                            <ColorMapEntry color="#007986" label="6.3" opacity="1.0" quantity="6.31579"/>
                            <ColorMapEntry color="#00a15e" label="8.4" opacity="1.0" quantity="8.42105"/>
                            <ColorMapEntry color="#00ca35" label="10.5" opacity="1.0" quantity="10.5263"/>
                            <ColorMapEntry color="#00f20d" label="12.6" opacity="1.0" quantity="12.6316"/>
                            <ColorMapEntry color="#1aff00" label="14.7" opacity="1.0" quantity="14.7368"/>
                            <ColorMapEntry color="#43ff00" label="16.8" opacity="1.0" quantity="16.8421"/>
                            <ColorMapEntry color="#6bff00" label="18.9" opacity="1.0" quantity="18.9474"/>
                            <ColorMapEntry color="#94ff00" label="21.1" opacity="1.0" quantity="21.0526"/>
                            <ColorMapEntry color="#bcff00" label="23.2" opacity="1.0" quantity="23.1579"/>
                            <ColorMapEntry color="#e5ff00" label="25.3" opacity="1.0" quantity="25.2632"/>
                            <ColorMapEntry color="#fff200" label="27.4" opacity="1.0" quantity="27.3684"/>
                            <ColorMapEntry color="#ffca00" label="29.5" opacity="1.0" quantity="29.4737"/>
                            <ColorMapEntry color="#ffa100" label="31.6" opacity="1.0" quantity="31.5789"/>
                            <ColorMapEntry color="#ff7900" label="33.7" opacity="1.0" quantity="33.6842"/>
                            <ColorMapEntry color="#ff5000" label="35.8" opacity="1.0" quantity="35.7895"/>
                            <ColorMapEntry color="#ff2800" label="37.9" opacity="1.0" quantity="37.9"/>
                            <ColorMapEntry color="#ff0000" label="100.0" opacity="1.0" quantity="100"/>
                            <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
                        </ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
