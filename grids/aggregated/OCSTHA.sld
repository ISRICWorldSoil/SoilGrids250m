<?xml version="1.0" ?>
<StyledLayerDescriptor version="1.0.0" xmlns="http://www.opengis.net/sld" xmlns:gml="http://www.opengis.net/gml" xmlns:ogc="http://www.opengis.net/ogc" xmlns:sld="http://www.opengis.net/sld">
    <UserLayer>
        <LayerFeatureConstraints>
            <FeatureTypeConstraint/>
        </LayerFeatureConstraints>
        <UserStyle>
            <Name>OCSTHA</Name>
            <Title/>
            <FeatureTypeStyle>
                <Name/>
                <Rule>
                    <RasterSymbolizer>
                     <Geometry>
                      <PropertyName>Soil organic carbon stock in tonnes per ha</PropertyName>
                      <Opacity>1</Opacity>
                    </Geometry>
                        <ColorMap>
                            <ColorMapEntry color="#FFFFFF" quantity="-32768" label="NODATA" opacity="0.0"/>
                            <ColorMapEntry color="#fff7f3" label="0" opacity="1.0" quantity="0"/>
                            <ColorMapEntry color="#fde0dd" label="56" opacity="1.0" quantity="56.25"/>
                            <ColorMapEntry color="#fcc5c0" label="113" opacity="1.0" quantity="112.5"/>
                            <ColorMapEntry color="#fa9fb5" label="169" opacity="1.0" quantity="168.75"/>
                            <ColorMapEntry color="#f768a1" label="225" opacity="1.0" quantity="225"/>
                            <ColorMapEntry color="#dd3497" label="281" opacity="1.0" quantity="281.25"/>
                            <ColorMapEntry color="#ae017e" label="338" opacity="1.0" quantity="337.5"/>
                            <ColorMapEntry color="#7a0177" label="394" opacity="1.0" quantity="393.75"/>
                            <ColorMapEntry color="#49006a" label="450" opacity="1.0" quantity="450"/>
                        </ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
