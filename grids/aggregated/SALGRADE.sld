<?xml version="1.0" ?>
<StyledLayerDescriptor version="1.0.0" xmlns="http://www.opengis.net/sld" xmlns:gml="http://www.opengis.net/gml" xmlns:ogc="http://www.opengis.net/ogc" xmlns:sld="http://www.opengis.net/sld">
    <UserLayer>
        <LayerFeatureConstraints>
            <FeatureTypeConstraint/>
        </LayerFeatureConstraints>
        <UserStyle>
            <Name>SLGWRB</Name>
            <Title/>
            <FeatureTypeStyle>
                <Name/>
                <Rule>
                    <RasterSymbolizer>
                     <Geometry>
                      <PropertyName>Saline / sodic soil grade based on WRB soil types and soil pH (SoilGrids250m)</PropertyName>
                      <Opacity>1</Opacity>
                    </Geometry>
                        <sld:ColorMap>
                            <sld:ColorMapEntry color="#f7f4f9" label="No saline/sodic properties" opacity="1.0" quantity="0"/>
                            <sld:ColorMapEntry color="#d4b9da" label="Some sodic properties" opacity="1.0" quantity="1"/>
                            <sld:ColorMapEntry color="#df65b0" label="Moderate sodic and solodic properties" opacity="1.0" quantity="2"/>
                            <sld:ColorMapEntry color="#ce1256" label="Solonetz or Calcisols, sodic properties and/or ph &gt;8.1" opacity="1.0" quantity="3"/>
                            <sld:ColorMapEntry color="#67001f" label="Solonetz, Solonchak, gypsic properties and/or ph &gt;8.5" opacity="1.0" quantity="4"/>
                        </sld:ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
