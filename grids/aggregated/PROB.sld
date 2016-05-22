<?xml version="1.0" ?>
<StyledLayerDescriptor version="1.0.0" xmlns="http://www.opengis.net/sld" xmlns:gml="http://www.opengis.net/gml" xmlns:ogc="http://www.opengis.net/ogc" xmlns:sld="http://www.opengis.net/sld">
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
                            <ogc:PropertyName>grid</ogc:PropertyName>
                        </Geometry>
                        <Opacity>1</Opacity>
                        <ColorMap>
                            <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
                            <ColorMapEntry color="#ffffd4" label="0.0" opacity="1.0" quantity="0"/>
                            <ColorMapEntry color="#fef7c5" label="2.1" opacity="1.0" quantity="2.10526"/>
                            <ColorMapEntry color="#feefb6" label="4.2" opacity="1.0" quantity="4.21053"/>
                            <ColorMapEntry color="#fee7a7" label="6.3" opacity="1.0" quantity="6.31579"/>
                            <ColorMapEntry color="#fedf99" label="8.4" opacity="1.0" quantity="8.42105"/>
                            <ColorMapEntry color="#fed588" label="10.5" opacity="1.0" quantity="10.5263"/>
                            <ColorMapEntry color="#fec873" label="12.6" opacity="1.0" quantity="12.6316"/>
                            <ColorMapEntry color="#feba5e" label="14.7" opacity="1.0" quantity="14.7368"/>
                            <ColorMapEntry color="#fead48" label="16.8" opacity="1.0" quantity="16.8421"/>
                            <ColorMapEntry color="#fe9f33" label="18.9" opacity="1.0" quantity="18.9474"/>
                            <ColorMapEntry color="#fa9226" label="21.1" opacity="1.0" quantity="21.0526"/>
                            <ColorMapEntry color="#f28620" label="23.2" opacity="1.0" quantity="23.1579"/>
                            <ColorMapEntry color="#ea7a1a" label="25.3" opacity="1.0" quantity="25.2632"/>
                            <ColorMapEntry color="#e26e15" label="27.4" opacity="1.0" quantity="27.3684"/>
                            <ColorMapEntry color="#da620f" label="29.5" opacity="1.0" quantity="29.4737"/>
                            <ColorMapEntry color="#ce580c" label="31.6" opacity="1.0" quantity="31.5789"/>
                            <ColorMapEntry color="#c14f0a" label="33.7" opacity="1.0" quantity="33.6842"/>
                            <ColorMapEntry color="#b34608" label="35.8" opacity="1.0" quantity="35.7895"/>
                            <ColorMapEntry color="#a63d06" label="37.9" opacity="1.0" quantity="37.8947"/>
                            <ColorMapEntry color="#993404" label="40.0" opacity="1.0" quantity="40"/>
                        </ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
